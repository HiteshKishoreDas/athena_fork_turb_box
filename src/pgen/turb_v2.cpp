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
static Real total_cooling;

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

static int cooling_flag = 0;

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
bool cloud_init_flag   = true;


void read_input (ParameterInput *pin){
  /*
  Reads all the static variables from the input file
  */

  n1 = pin->GetInteger("mesh","nx1");
  n2 = pin->GetInteger("mesh","nx2");
  n3 = pin->GetInteger("mesh","nx3");

  cooling_flag = pin->GetInteger("problem","cooling_flag");

  amb_rho      = pin->GetReal("problem","amb_rho");

  cloud_radius = pin->GetReal("problem","cloud_radius");
  cloud_time   = pin->GetReal("problem","cloud_time");
  cloud_chi    = pin->GetReal("problem","cloud_chi");

  heating_rate = pin->GetReal("problem","heating");

  T_floor      = pin->GetReal("problem","T_floor");
  T_ceil       = pin->GetReal("problem","T_ceil");
  T_hot        = pin->GetReal("problem","T_hot");
  T_hot_req    = pin->GetReal("problem","T_hot_req");
  T_cold       = pin->GetReal("problem","T_cold");
  T_cut_mul    = pin->GetReal("problem","T_cut_mul");
  T_cut        = pin->GetReal("problem","T_cut");

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

              cons(IEN,k,j,i) += ccool + heating_rate*dt;
              total_cooling -= ccool;
            }


            //*_____________________________
            


            // ________________________________
            // FOR DEBUG PURPOSES
            // if ((temp>T_cut) && ((temp_new-temp)!=0.0)){

                // printf("T, delT: %lf %.60lf\n", temp, temp_new-temp);
            // }
            // ________________________________


          
        }

        else{  // If T<=T_floor

          printf("+");

          Real ccool = ((T_floor-temp)/(KELVIN*mu))*cons(IDN,k,j,i)/(g-1); 
          cons(IEN,k,j,i) += ccool;
          total_cooling -= ccool;

        }

      }
    }
  }

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

  if (cooling_flag!=0){
  
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

    AllocateUserHistoryOutput(11);

    EnrollUserHistoryOutput(0, rho_sum, "rho_sum");
    EnrollUserHistoryOutput(1, rho_sq_sum, "rho_sq_sum");
    EnrollUserHistoryOutput(2, c_s_sum, "c_s_sum");
    EnrollUserHistoryOutput(3, Pth_sum, "Pth_sum");
    EnrollUserHistoryOutput(4, PB_sum, "PB_sum");
    EnrollUserHistoryOutput(5, Bx_sum, "Bx_sum");
    EnrollUserHistoryOutput(6, By_sum, "By_sum");
    EnrollUserHistoryOutput(7, Bz_sum, "Bz_sum");
    EnrollUserHistoryOutput(8, cold_gas, "cold_gas");
    EnrollUserHistoryOutput(9, hst_total_cooling, "total_cooling");
    EnrollUserHistoryOutput(10, T_sum, "T_sum");

  }
  else {
    AllocateUserHistoryOutput(6);

    EnrollUserHistoryOutput(0, rho_sum, "rho_sum");
    EnrollUserHistoryOutput(1, rho_sq_sum, "rho_sq_sum");
    EnrollUserHistoryOutput(2, c_s_sum, "c_s_sum");
    EnrollUserHistoryOutput(3, cold_gas, "cold_gas");
    EnrollUserHistoryOutput(4, hst_total_cooling, "total_cooling");
    EnrollUserHistoryOutput(5, T_sum, "T_sum");
  }

  
  return;
}




//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  Real g   = pin->GetReal("hydro","gamma");

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

  return;

}


void MeshBlock::InitUserMeshBlockData(ParameterInput *pin){
  AllocateUserOutputVariables(1);
  return;
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin){

  Real g = peos->GetGamma();

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
      
      }
    }
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
