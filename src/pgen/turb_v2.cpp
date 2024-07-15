//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file turb.cpp
//! \brief Problem generator for turbulence driver

// C headers

// C++ headers
#include <cmath>
#include <ctime>
#include <sstream>
#include <stdexcept>

#include<memory>    // std::unique_ptr

//! Remove if not required
#include <algorithm>     // max()
#include <string>     // c_str(), string

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../fft/athena_fft.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../utils/utils.hpp"
#include "../scalars/scalars.hpp"
#include "../inputs/hdf5_reader.hpp"  // HDF5ReadRealArray()

//! Remove if not required
// #include "../inputs/hdf5_reader.hpp"  // HDF5ReadRealArray()

// #include "../utils/townsend_cooling.hpp"  // T_new() : Implementation by Hitesh (buggy)
#include "../utils/townsend_cooling_max.hpp" // T_new() : Impolementation by Max (works)

#include "../utils/hst_func.hpp"             // All history output functions
#include "../utils/code_units.hpp"           // Code units and constants

//* ___________________________
//* For Max's townsend cooling

using namespace std;

// Helper functions
void myprint(string msg);


static Real tfloor, tnotcool, tcut_hst, r_drop;
static Real Lambda_fac, Lambda_fac_time;         // for boosting cooling
static Real total_cooling = 0.0;

// Returns unique pointer
// This is a function in C++13 onwards
// The implementation here is copied from https://stackoverflow.com/a/17903225/1834164

template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args)
{
  return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

std::unique_ptr<Cooling> cooler;


//*_____________________________


#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

// User headers
#include <fstream>    // for file io
#include <iomanip>    // for std::setprecision


// Variable read from the input file

static int n1 = 64;
static int n2 = 64;
static int n3 = 64;

static Real pfloor = 0;

static int cooling_flag = 0;
static int global_cooling_flag = 0;

static Real amb_rho = 1.0;

static Real cloud_radius;
static Real cloud_time;
static Real cloud_chi;

static Real heating_rate;

//TODO: Double definition at code_unit.hpp, figure out a better way

// Cloud position
static Real cloud_pos_x = 0.0;
static Real cloud_pos_y = 0.0;
static Real cloud_pos_z = 0.0;

// Magnetic fields
static Real B_x = 0.0;
static Real B_y = 0.0;
static Real B_z = 1.0;

// static Real mH = 1.0; //! Not required?

static bool cooling_flag_print_count = false;

bool temp_rescale_flag = true;
int hdf_cloud_flag   = 0;
int rescale_flag = 0;

static std::string b1_input_filename;
static std::string b2_input_filename;
static std::string b3_input_filename;

static std::string dataset_b1;
static std::string dataset_b2;
static std::string dataset_b3;

void read_input (ParameterInput *pin){
  /*
  Reads all the static variables from the input file
  */

  n1 = pin->GetInteger("mesh","nx1");
  n2 = pin->GetInteger("mesh","nx2");
  n3 = pin->GetInteger("mesh","nx3");

  cooling_flag = pin->GetInteger("problem","cooling_flag");
  global_cooling_flag = pin->GetInteger("problem","global_cooling_flag");
  hdf_cloud_flag = pin->GetInteger("problem","hdf_cloud_flag");
  rescale_flag = pin->GetInteger("problem","rescale_flag");

  amb_rho      = pin->GetReal("problem","amb_rho");

  cloud_radius = pin->GetReal("problem","cloud_radius");
  cloud_time   = pin->GetReal("problem","cloud_time");
  cloud_chi    = pin->GetReal("problem","cloud_chi");

  heating_rate = pin->GetReal("problem","heating");

  pfloor       = pin->GetReal("hydro","pfloor"); // get the pressure floor
  T_floor      = pin->GetReal("problem","T_floor");
  T_ceil       = pin->GetReal("problem","T_ceil");
  T_hot        = pin->GetReal("problem","T_hot");
  T_hot_req    = pin->GetReal("problem","T_hot_req");
  T_cold       = pin->GetReal("problem","T_cold");
  T_warm       = pin->GetReal("problem","T_warm");
  T_cut_mul    = pin->GetReal("problem","T_cut_mul");
  T_cut        = pin->GetReal("problem","T_cut");

  printf("pfloor = %lf\n", pfloor);
  printf("T_floor = %lf\n", T_floor);
  printf("T_ceil = %lf\n", T_ceil);
  printf("T_hot = %lf\n", T_hot);

  // T_cut = T_cut_mul*T_hot_req;

  // Xsol            = pin->GetReal("problem","Xsol");
  // Zsol            = pin->GetReal("problem","Zsol");

  // X = Xsol * 0.7381;
  // Z = Zsol * 0.0134;
  // Y = 1 - X - Z;

  // mu  = 1.0/(2.*X+ 3.*(1.-X-Z)/4.+ Z/2.);
  // mue = 2.0/(1.0+X);
  // muH = 1.0/X;

  cloud_pos_x  = pin->GetReal("problem","cloud_pos_x");
  cloud_pos_y  = pin->GetReal("problem","cloud_pos_y");
  cloud_pos_z  = pin->GetReal("problem","cloud_pos_z");

  if (MAGNETIC_FIELDS_ENABLED) {
    B_x = pin->GetReal("problem", "B_x");
    B_y = pin->GetReal("problem", "B_y");
    B_z = pin->GetReal("problem", "B_z");
  }

  printf("____ Input file read! ______\n");

  return;

}


void townsend_cooling(MeshBlock *pmb, const Real time, const Real dt,
             const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
             AthenaArray<Real> &cons_scalar) {

  Real g = pmb->peos->GetGamma();

  Real sim_time = pmb->pmy_mesh->time;

  Real t_cloud = cloud_time;

  // printf("T_floor =  %lf\n", T_floor);

  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    for (int j = pmb->js; j <= pmb->je; ++j) {
      for (int i = pmb->is; i <= pmb->ie; ++i) {

        Real temp = (prim(IPR,k,j,i) / cons(IDN,k,j,i)) * KELVIN * mu ;

        if (temp > T_floor) {

          if (cooling_flag!=0){

            //* For own townsend cooling

            // Real *lam_para = Lam_file_read(temp);

            //printf("%f %f %f\n",lam_para[0],lam_para[1],lam_para[2]);

            // Real temp_new = T_new(temp, lam_para[0], lam_para[1], lam_para[2],
            //                         prim(IDN,k,j,i), dt, 
            //                         mu, mue, muH, 
            //                         g, T_floor, T_ceil, T_cut);

                            // Defined in ../utils/townsend_cooling.hpp

            //  NOTE: Above, dt is in code units to avoid overflow. unit_time is cancelled in 
            //  calculation of T_new as we calculate (dt/tcool)

            // cons(IEN,k,j,i) += ((temp_new-temp)/(KELVIN*mu))*cons(IDN,k,j,i)/(g-1);
            
            //* For Max's townsend cooling

            if (temp<T_cut){
              Real rho      = cons(IDN,k,j,i);
              Real temp_cgs = temp;
              Real rho_cgs  = rho * unit_density;
              Real dt_cgs   = dt  * unit_time;
              Real cLfac    = 1.0;

              Real temp_new = max(cooler->townsend(temp_cgs,rho_cgs,dt_cgs, 1.0), T_floor);

              Real ccool = ((temp_new-temp)/(KELVIN*mu))*cons(IDN,k,j,i)/(g-1);

              
              cons(IEN,k,j,i) += ccool;
      

              // printf("temp = %lf \n", temp);
              // printf("temp_new = %lf \n", temp_new);
              total_cooling -= ccool;
            } // if (temp<T_cut)

            // ________________________________
            // FOR DEBUG PURPOSES
            // if ((temp>T_cut) && ((temp_new-temp)!=0.0)){

                // printf("T, delT: %lf %.60lf\n", temp, temp_new-temp);
            // }
            // ________________________________


          } // if (cooling_flag!=0)
          
        } // if (temp > T_floor)

        else{  // If T<=T_floor

          // printf("+");

          Real ccool = ((T_floor-temp)/(KELVIN*mu))*cons(IDN,k,j,i)/(g-1); 
          cons(IEN,k,j,i) += ccool;
          // total_cooling -= ccool;

        } // else, If T<=T_floor

        if (global_cooling_flag!=0){
         cons(IEN,k,j,i) += heating_rate*dt;
        }

      } // for (int i = pmb->is; i <= pmb->ie; ++i)
    } // for (int j = pmb->js; j <= pmb->je; ++j)
  } // for (int k = pmb->ks; k <= pmb->ke; ++k)

  return;
}


Real hst_total_cooling(MeshBlock *pmb, int iout) {
  if(pmb->lid == 0)
    return total_cooling;
  else
    return 0;
}

//* Source Terms
void Source(MeshBlock *pmb, const Real time, const Real dt,
             const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
             AthenaArray<Real> &cons_scalar) {

  if ((cooling_flag!=0)||(global_cooling_flag!=0)){
  
    if (!cooling_flag_print_count){ 
        printf("___________________________________________\n");
        printf("!! Cooling included .......................\n");
        printf("___________________________________________\n");
    
        cooling_flag_print_count = true;
    }
  
    townsend_cooling(pmb, time, dt,
            prim, prim_scalar,
            bcc, cons,
            cons_scalar);
  }

  return;

}




//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief
//========================================================================================
void Mesh::InitUserMeshData(ParameterInput *pin) {


  if (SELF_GRAVITY_ENABLED) {
    Real four_pi_G = pin->GetReal("problem","four_pi_G");
    Real eps = pin->GetOrAddReal("problem","grav_eps", 0.0);
    SetFourPiG(four_pi_G);
    SetGravityThreshold(eps);
  }

  // turb_flag is initialzed in the Mesh constructor to 0 by default;
  // turb_flag = 1 for decaying turbulence
  // turb_flag = 2 for impulsively driven turbulence
  // turb_flag = 3 for continuously driven turbulence
  turb_flag = pin->GetInteger("problem","turb_flag");
  if (turb_flag != 0) {
    #ifndef FFT
        std::stringstream msg;
        msg << "### FATAL ERROR in TurbulenceDriver::TurbulenceDriver" << std::endl
            << "non zero Turbulence flag is set without FFT!" << std::endl;
        ATHENA_ERROR(msg);
        return;
    #endif
  }

  //* Read static variables from the input file
  read_input(pin);

  //* Enroll the Source terms
  EnrollUserExplicitSourceFunction(Source);
  
  cooler = make_unique<Cooling>();
  // std::unique_ptr<Cooling> cooler = std::make_unique<Cooling>();
  // cooler = std::make_unique(Cooling);

  //* History outputs
  if (MAGNETIC_FIELDS_ENABLED) {

    AllocateUserHistoryOutput(12);

    EnrollUserHistoryOutput(0, rho_sum, "rho_sum");
    EnrollUserHistoryOutput(1, rho_sq_sum, "rho_sq_sum");
    EnrollUserHistoryOutput(2, c_s_sum, "c_s_sum");
    EnrollUserHistoryOutput(3, Pth_sum, "Pth_sum");
    EnrollUserHistoryOutput(4, PB_sum, "PB_sum");
    EnrollUserHistoryOutput(5, Bx_sum, "Bx_sum");
    EnrollUserHistoryOutput(6, By_sum, "By_sum");
    EnrollUserHistoryOutput(7, Bz_sum, "Bz_sum");
    EnrollUserHistoryOutput(8, cold_gas, "cold_gas");
    EnrollUserHistoryOutput(9, warm_gas, "warm_gas");
    EnrollUserHistoryOutput(10, hst_total_cooling, "total_cooling");
    EnrollUserHistoryOutput(11, T_sum, "T_sum");

  }
  else {
    AllocateUserHistoryOutput(7);

    EnrollUserHistoryOutput(0, rho_sum, "rho_sum");
    EnrollUserHistoryOutput(1, rho_sq_sum, "rho_sq_sum");
    EnrollUserHistoryOutput(2, c_s_sum, "c_s_sum");
    EnrollUserHistoryOutput(3, cold_gas, "cold_gas");
    EnrollUserHistoryOutput(4, warm_gas, "warm_gas");
    EnrollUserHistoryOutput(5, hst_total_cooling, "total_cooling");
    EnrollUserHistoryOutput(6, T_sum, "T_sum");
  }

  
  return;
}




//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  Real g   = pin->GetReal("hydro","gamma");
  constexpr Real scalar_norm = NSCALARS > 0 ? (Real)NSCALARS : 1.0;

  // Real l_shatter = cooler->tcool(tfloor * unit_temp,
	// 			 drat * t_drop / tfloor * unit_rho) / unit_time * \
  //   sqrt(gamma * tfloor);

  // Read/set problem parameters
  Real sound_speed = sqrt(g*CONST_kB*T_floor / (mu*CONST_amu) );
  Real l_shatter = cooler->tcool(  T_floor, cloud_chi* amb_rho * unit_density) ;
  l_shatter     *= sound_speed ;
  l_shatter     /= unit_length;


  //! Change t_cool,cold to t_cool,floor
  //! T_floor is different from T_cold (=2*T_floor usually), which is a cutoff temp
  //! Can be very confusing...
  stringstream msg;
  if(gid == 0) {
    msg << "[pgen] unit_time = " << unit_time << ", unit_rho = " << unit_density 
        << ", unit_len = " << unit_length;
    myprint(msg.str());
    msg.str("");
    msg << "[pgen] T_floor = " << T_floor << ", T_cut = " << T_cut << ", T_hot = " << T_hot
        << ", t_cool,hot = " << cooler->tcool(T_hot,
                                              amb_rho * unit_density) / unit_time
        << ", t_cool,mix = " << cooler->tcool(sqrt(T_hot*T_cold),
                                              sqrt(cloud_chi)* amb_rho * unit_density) / unit_time
        << ", t_cool,cold = " << cooler->tcool(T_cold,
                                              cloud_chi* amb_rho * unit_density) / unit_time
        << ", sound_speed = " << sound_speed 
        << ", R_cloud = " <<  cloud_radius
        << ", l_shatter = " << l_shatter;


    myprint(msg.str());
    msg.str("");
  } 


  // * Read static variables from the input file
  read_input(pin);

  if (hdf_cloud_flag==0){
    printf("Before pgen loop...................\n");

    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          phydro->u(IDN,k,j,i) = amb_rho;

          phydro->u(IM1,k,j,i) = 0.0;
          phydro->u(IM2,k,j,i) = 0.0;
          phydro->u(IM3,k,j,i) = 0.0;

          if (NON_BAROTROPIC_EOS) {
            phydro->u(IEN,k,j,i) = (T_hot/(KELVIN*mu))*phydro->u(IDN,k,j,i)/(g-1);
          }

          // Add scalar to the cloud
          if (NSCALARS > 0) {
            for (int n=0; n<NSCALARS; ++n) {
              pscalars->s(n,k,j,i)  = 0.0 * phydro->u(IDN, k, j, i)/scalar_norm;
          }
          }
          // printf("rho: %f\n",phydro->u(IDN,k,j,i));
        }
      }
    }

    printf("After pgen loop...................\n");

    // initialize interface B
    if (MAGNETIC_FIELDS_ENABLED) {

      printf("________________________________________________\n");
      printf("Inside Magnetic initialisation \n");
      printf("________________________________________________\n");

      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          for (int i=is; i<=ie+1; ++i) {
            pfield->b.x1f(k,j,i) = B_x;
          }
        }
      }
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je+1; ++j) {
          for (int i=is; i<=ie; ++i) {
            pfield->b.x2f(k,j,i) = B_y;
          }
        }
      }
      for (int k=ks; k<=ke+1; ++k) {
        for (int j=js; j<=je; ++j) {
          for (int i=is; i<=ie; ++i) {
            pfield->b.x3f(k,j,i) = B_z;
          }
        }
      }
      if (NON_BAROTROPIC_EOS) {
        for (int k=ks; k<=ke; ++k) {
          for (int j=js; j<=je; ++j) {
            for (int i=is; i<=ie; ++i) {
              phydro->u(IEN,k,j,i) += 0.5*(SQR(B_x) + SQR(B_y) + SQR(B_z));
            }
          }
        }
      }
    }

    printf("End of pgen.............................\n");
  }
  else { // hdf_cloud_flag == 1

    // Determine locations of initial values
    std::string input_filename = pin->GetString("problem", "cons_input_filename");
    std::string dataset_cons = pin->GetString("problem", "dataset_cons");
    int index_dens = pin->GetInteger("problem", "index_dens");
    int index_mom1 = pin->GetInteger("problem", "index_mom1");
    int index_mom2 = pin->GetInteger("problem", "index_mom2");
    int index_mom3 = pin->GetInteger("problem", "index_mom3");
    int index_etot = pin->GetInteger("problem", "index_etot");

  
    // Set conserved array selections
    int start_cons_file[5];
    start_cons_file[1] = gid;
    start_cons_file[2] = 0;
    start_cons_file[3] = 0;
    start_cons_file[4] = 0;
    int start_cons_indices[5];
    start_cons_indices[IDN] = index_dens;
    start_cons_indices[IM1] = index_mom1;
    start_cons_indices[IM2] = index_mom2;
    start_cons_indices[IM3] = index_mom3;
    start_cons_indices[IEN] = index_etot;
    int count_cons_file[5];
    count_cons_file[0] = 1;
    count_cons_file[1] = 1;
    count_cons_file[2] = block_size.nx3;
    count_cons_file[3] = block_size.nx2;
    count_cons_file[4] = block_size.nx1;
    int start_cons_mem[4];
    start_cons_mem[1] = ks;
    start_cons_mem[2] = js;
    start_cons_mem[3] = is;
    int count_cons_mem[4];
    count_cons_mem[0] = 1;
    count_cons_mem[1] = block_size.nx3;
    count_cons_mem[2] = block_size.nx2;
    count_cons_mem[3] = block_size.nx1;
  
    // Set conserved values from file
    for (int n = 0; n < NHYDRO; ++n) {
      start_cons_file[0] = start_cons_indices[n];
      start_cons_mem[0] = n;
      HDF5ReadRealArray(input_filename.c_str(), dataset_cons.c_str(), 5, start_cons_file,
                        count_cons_file, 4, start_cons_mem,
                        count_cons_mem, phydro->u, true);
    }

    // Set field array selections
    int start_field_file[4];
    start_field_file[0] = gid;
    start_field_file[1] = 0;
    start_field_file[2] = 0;
    start_field_file[3] = 0;
    int count_field_file[4];
    count_field_file[0] = 1;
    int start_field_mem[4];
    start_field_mem[0] = 0;
    start_field_mem[1] = ks;
    start_field_mem[2] = js;
    start_field_mem[3] = is;
    int count_field_mem[3];
  
  
    // Set magnetic field values from file
    if (MAGNETIC_FIELDS_ENABLED) {
      b1_input_filename = pin->GetString("problem", "b1_input_filename");
      b2_input_filename = pin->GetString("problem", "b2_input_filename");
      b3_input_filename = pin->GetString("problem", "b3_input_filename");

      dataset_b1 = "hydro";//pin->GetString("problem", "dataset_b1");
      dataset_b2 = "hydro";//pin->GetString("problem", "dataset_b2");
      dataset_b3 = "hydro";//pin->GetString("problem", "dataset_b3");


      // Set B1
      count_field_file[1] = block_size.nx3;
      count_field_file[2] = block_size.nx2;
      count_field_file[3] = block_size.nx1 + 1;
      count_field_mem[0] = block_size.nx3;
      count_field_mem[1] = block_size.nx2;
      count_field_mem[2] = block_size.nx1 + 1;
      HDF5ReadRealArray(b1_input_filename.c_str(), dataset_b1.c_str(), 5, start_field_file,
                        count_field_file, 4, start_field_mem,
                        count_field_mem, pfield->b.x1f, true);
  
      // Set B2
      count_field_file[1] = block_size.nx3;
      count_field_file[2] = block_size.nx2 + 1;
      count_field_file[3] = block_size.nx1;
      count_field_mem[0] = block_size.nx3;
      count_field_mem[1] = block_size.nx2 + 1;
      count_field_mem[2] = block_size.nx1;
      HDF5ReadRealArray(b2_input_filename.c_str(), dataset_b2.c_str(), 5, start_field_file,
                        count_field_file, 4, start_field_mem,
                        count_field_mem, pfield->b.x2f, true);
  
      // Set B3
      count_field_file[1] = block_size.nx3 + 1;
      count_field_file[2] = block_size.nx2;
      count_field_file[3] = block_size.nx1;
      count_field_mem[0] = block_size.nx3 + 1;
      count_field_mem[1] = block_size.nx2;
      count_field_mem[2] = block_size.nx1;
      HDF5ReadRealArray(b3_input_filename.c_str(), dataset_b3.c_str(), 5, start_field_file,
                        count_field_file, 4, start_field_mem,
                        count_field_mem, pfield->b.x3f, true);
    }
  
    // Make no-op collective reads if using MPI and ranks have unequal numbers of blocks
  #ifdef MPI_PARALLEL
    {
      int num_blocks_this_rank = pmy_mesh->nblist[Globals::my_rank];
      if (lid == num_blocks_this_rank - 1) {
        int block_shortage_this_rank = 0;
        for (int rank = 0; rank < Globals::nranks; ++rank) {
          block_shortage_this_rank =
              std::max(block_shortage_this_rank,
                      pmy_mesh->nblist[rank] - num_blocks_this_rank);
        }
        for (int block = 0; block < block_shortage_this_rank; ++block) {
          for (int n = 0; n < NHYDRO; ++n) {
            start_cons_file[0] = start_cons_indices[n];
            start_cons_mem[0] = n;
            HDF5ReadRealArray(input_filename.c_str(), dataset_cons.c_str(), 5,
                              start_cons_file, count_cons_file, 4,
                              start_cons_mem, count_cons_mem,
                              phydro->u, true, true);
          }
          if (MAGNETIC_FIELDS_ENABLED) {
            count_field_file[1] = block_size.nx3;
            count_field_file[2] = block_size.nx2;
            count_field_file[3] = block_size.nx1 + 1;
            count_field_mem[0] = block_size.nx3;
            count_field_mem[1] = block_size.nx2;
            count_field_mem[2] = block_size.nx1 + 1;
            HDF5ReadRealArray(b1_input_filename.c_str(), dataset_b1.c_str(), 5,
                              start_field_file, count_field_file, 4,
                              start_field_mem, count_field_mem,
                              pfield->b.x1f, true, true);
            count_field_file[1] = block_size.nx3;
            count_field_file[2] = block_size.nx2 + 1;
            count_field_file[3] = block_size.nx1;
            count_field_mem[0] = block_size.nx3;
            count_field_mem[1] = block_size.nx2 + 1;
            count_field_mem[2] = block_size.nx1;
            HDF5ReadRealArray(b2_input_filename.c_str(), dataset_b2.c_str(), 5,
                              start_field_file, count_field_file, 4,
                              start_field_mem, count_field_mem,
                              pfield->b.x2f, true, true);
            count_field_file[1] = block_size.nx3 + 1;
            count_field_file[2] = block_size.nx2;
            count_field_file[3] = block_size.nx1;
            count_field_mem[0] = block_size.nx3 + 1;
            count_field_mem[1] = block_size.nx2;
            count_field_mem[2] = block_size.nx1;
            HDF5ReadRealArray(b3_input_filename.c_str(), dataset_b3.c_str(), 5,
                              start_field_file, count_field_file, 4,
                              start_field_mem, count_field_mem,
                              pfield->b.x3f, true, true);
          }
        }
      }
    }
  #endif


  //* Rescale temperatures
  if (rescale_flag==1){

    Real local_T_sum_1 = 0.0;
    Real global_T_sum_1;

    // Calculate average temperature before rescaling
    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je; ++j) {
        for (int i = is; i <= ie; ++i) {
          Real KE = 0.0;
          KE += phydro->u(IM1,k,j,i) * phydro->u(IM1,k,j,i);
          KE += phydro->u(IM2,k,j,i) * phydro->u(IM2,k,j,i);
          KE += phydro->u(IM3,k,j,i) * phydro->u(IM3,k,j,i);
          KE /= 2.0*phydro->u(IDN,k,j,i);

          Real prs = phydro->u(IEN,k,j,i);
          prs -= KE;

          if (MAGNETIC_FIELDS_ENABLED) {
            Real BE = 0.0;
            BE += 0.5 * pfield->b.x1f(k,j,i) * pfield->b.x1f(k,j,i);
            BE += 0.5 * pfield->b.x2f(k,j,i) * pfield->b.x2f(k,j,i);
            BE += 0.5 * pfield->b.x3f(k,j,i) * pfield->b.x3f(k,j,i);

            prs -= BE;
          }

          prs *= g-1;

          Real temp  = (prs / phydro->u(IDN,k,j,i)) * KELVIN * mu ;

          local_T_sum_1 += temp;

        }
      }
    }//End of for loop over domain

    #ifdef MPI_PARALLEL //Parallel
      MPI_Allreduce(&local_T_sum_1, &global_T_sum_1, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
      Real N_cells = pmy_mesh->nbtotal*(ke-ks+1)*(je-js+1)*(ie-is+1);
    #else  // Serial
      global_T_sum_1 = local_T_sum_1;
      Real N_cells = (ke-ks+1)*(je-js+1)*(ie-is+1);
    #endif // #ifdef MPI_PARALLEL

    Real T_avg = global_T_sum_1/N_cells; //* Average temperature
    printf("T_avg: %lf \n",T_avg);
    printf("local_T_sum_1: %lf \n",local_T_sum_1);
    printf("local avg T: %lf \n",local_T_sum_1 / 32768.);  // individual cell temperature
    printf("global_T_sum_1: %lf \n",global_T_sum_1);
    printf("N_cells: %lf \n",N_cells);
    printf("scaling_factor: %lf \n",T_hot_req/T_avg);

    // Rescaling the temperature
    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je; ++j) {
        for (int i = is; i <= ie; ++i) {
          Real KE = 0.0;
          KE += phydro->u(IM1,k,j,i) * phydro->u(IM1,k,j,i);
          KE += phydro->u(IM2,k,j,i) * phydro->u(IM2,k,j,i);
          KE += phydro->u(IM3,k,j,i) * phydro->u(IM3,k,j,i);
          KE /= 2.0*phydro->u(IDN,k,j,i);

          Real prs = phydro->u(IEN,k,j,i);
          prs -= KE;

          if (MAGNETIC_FIELDS_ENABLED) {
            Real BE = 0.0;
            BE += 0.5 * pfield->b.x1f(k,j,i) * pfield->b.x1f(k,j,i);
            BE += 0.5 * pfield->b.x2f(k,j,i) * pfield->b.x2f(k,j,i);
            BE += 0.5 * pfield->b.x3f(k,j,i) * pfield->b.x3f(k,j,i);

            prs -= BE;
          }

          prs *= g-1;

          Real temp  = (prs / phydro->u(IDN,k,j,i)) * KELVIN * mu ;
          Real T_new = temp*(T_hot_req/T_avg);

          phydro->u(IEN,k,j,i) += ((T_new-temp)/(KELVIN*mu))*phydro->u(IDN,k,j,i)/(g-1);

        }
      }
    } //End of for loop over domain


    printf("______________________________\n");
    printf("__Temperature rescaled!_______\n");
    printf("______________________________\n");

  } // End of temperature rescaling section


  for (int k = ks; k <= ke; ++k) {
    Real z_coord = pcoord->x3v(k);
    for (int j = js; j <= je; ++j) {
      Real y_coord = pcoord->x2v(j);
      for (int i = is; i <= ie; ++i) {
        Real x_coord = pcoord->x1v(i);

        Real dist = (x_coord-cloud_pos_x)*(x_coord-cloud_pos_x);
        dist     += (y_coord-cloud_pos_y)*(y_coord-cloud_pos_y);
        dist     += (z_coord-cloud_pos_z)*(z_coord-cloud_pos_z);

        // Inside the cloud region
        if(dist <= cloud_radius*cloud_radius){

          //Initial energy ______________________________________//
          Real TE_in = phydro->u(IEN,k,j,i);

          Real KE_in = phydro->u(IM1,k,j,i) * phydro->u(IM1,k,j,i);
          KE_in += phydro->u(IM2,k,j,i) * phydro->u(IM2,k,j,i);
          KE_in += phydro->u(IM3,k,j,i) * phydro->u(IM3,k,j,i);
          KE_in /= 2.0*phydro->u(IDN,k,j,i);

          Real BE_in = 0.0;
          if (MAGNETIC_FIELDS_ENABLED) {
            BE_in += 0.5 * pfield->b.x1f(k,j,i) * pfield->b.x1f(k,j,i);
            BE_in += 0.5 * pfield->b.x2f(k,j,i) * pfield->b.x2f(k,j,i);
            BE_in += 0.5 * pfield->b.x3f(k,j,i) * pfield->b.x3f(k,j,i);
          }

          Real IE_in = TE_in - KE_in - BE_in;

          //_______________________________________________________//

          Real rho_temp = phydro->u(IDN,k,j,i);

          Real temp = (IE_in * (g-1.0)/rho_temp) * KELVIN * mu ;

          Real chi_lim = std::min(cloud_chi,temp/T_floor);
          printf("chi_lim: %lf\n", chi_lim);


          phydro->u(IDN,k,j,i) *= chi_lim;//*2.0;

          phydro->u(IM1,k,j,i) /= sqrt(chi_lim);//*2.0;
          phydro->u(IM2,k,j,i) /= sqrt(chi_lim);//*2.0;
          phydro->u(IM3,k,j,i) /= sqrt(chi_lim);//*2.0;


          //Final energy _________________________________________//

          Real TE_fn = phydro->u(IEN,k,j,i);

          Real KE_fn = phydro->u(IM1,k,j,i) * phydro->u(IM1,k,j,i);
          KE_fn += phydro->u(IM2,k,j,i) * phydro->u(IM2,k,j,i);
          KE_fn += phydro->u(IM3,k,j,i) * phydro->u(IM3,k,j,i);
          KE_fn /= 2.0*phydro->u(IDN,k,j,i);

          Real BE_fn = 0.0;
          if (MAGNETIC_FIELDS_ENABLED) {
            BE_fn += 0.5 * pfield->b.x1f(k,j,i) * pfield->b.x1f(k,j,i);
            BE_fn += 0.5 * pfield->b.x2f(k,j,i) * pfield->b.x2f(k,j,i);
            BE_fn += 0.5 * pfield->b.x3f(k,j,i) * pfield->b.x3f(k,j,i);
          }

          Real IE_fn = TE_fn - KE_fn - BE_fn;

          //_______________________________________________________//

          phydro->u(IEN,k,j,i) = IE_in + KE_fn + BE_fn;


          // Add scalar to the cloud
          if (NSCALARS > 0) {
            for (int n=0; n<NSCALARS; ++n) {
              pscalars->s(n,k,j,i)  = 1.0*phydro->u(IDN,k,j,i)/scalar_norm;
            }
          }

          printf("Cloud added here! :) \n");
          printf("Density: %lf\n", phydro->u(IDN,k,j,i));

        } // If condition on cloud region closed

      } // End of loop over i
    } // End of loop over j
  } // End of loop over k

  printf("______________________________\n");
  printf("__Cloud added!________________\n");
  printf("______________________________\n");

  } 
    

  return;
}


void MeshBlock::InitUserMeshBlockData(ParameterInput *pin){
  AllocateUserOutputVariables(1);
  return;
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin){

  Real g = peos->GetGamma();

  // calcualte the number of cells for each mesh
  int num_cell_mesh = (ke - ks + 1) * (je - js + 1) * (ie - is + 1);
  int num_cell_pfloor = 0;  // number of cells < 2 * pfloor

  // loop through each cell in the mesh
  for(int k=ks; k<=ke; k++) {
    for(int j=js; j<=je; j++) {
      for(int i=is; i<=ie; i++) {

        Real lum_cell = 0;
        Real temp = (phydro->u(IPR,k,j,i) / phydro->u(IDN,k,j,i)) * KELVIN * mu;

        if (temp > T_floor) {

            //* For Max's townsend cooling

            if (temp<T_cut){
              Real rho      = phydro->u(IDN,k,j,i);
              Real temp_cgs = temp;
              Real rho_cgs  = rho * unit_density;
              Real dt_cgs   = pmy_mesh->dt * unit_time;
              Real cLfac    = 1.0;

              Real temp_new = max(cooler->townsend(temp_cgs,rho_cgs,dt_cgs, 1.0), T_floor);

              Real ccool = ((temp_new-temp)/(KELVIN*mu))*phydro->u(IDN,k,j,i)/(g-1);

              lum_cell-= ccool;
            }

            //*_____________________________
          
        }

        else{  // If T<=T_floor

          Real ccool = ((T_floor-temp)/(KELVIN*mu))*phydro->u(IDN,k,j,i)/(g-1); 
          lum_cell-= ccool;

        }
        
        user_out_var(0,k,j,i) = lum_cell/pmy_mesh->dt;      
      
        // evaluate if the cell hits pfloor
        if ((phydro->u(IPR,k,j,i)) < 2 * pfloor) {
          num_cell_pfloor++;
        }

      }
    }
    }

  // check for number of cells that are hitting the pressure floor
  if (num_cell_pfloor >= static_cast<float>(num_cell_mesh) * 0.005) {  // return an error if > 0.5%
    std::stringstream msg;
    msg << "### FATAL ERROR in function" << std::endl;
    std::cout << msg.str();
    std::exit(EXIT_FAILURE);
    return;
  }
}

//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//  \brief
//========================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
}


// prints message only once, and adds newline.
void myprint(string msg) {
  if (Globals::my_rank==0) {
    cout << msg << endl;
  }
}
