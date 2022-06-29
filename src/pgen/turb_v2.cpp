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

//! Remove if not required
// #include <algorithm>  // max()
// #include <string>     // c_str(), string

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

#include "../utils/townsend_cooling.hpp"  // T_new()
#include "../utils/hst_func.hpp"          // All history output functions

#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

// User headers
#include <fstream>    // for file io
#include <iomanip>    // for std::setprecision

// * Defined in townsend_cooling.hpp
//TODO: Maybe move the units to another file from townsend_cooling.hpp
/*
static const Real CONST_pc  = 3.086e18;
static const Real CONST_yr  = 3.154e7;
static const Real CONST_amu = 1.66053886e-24;
static const Real CONST_kB  = 1.3806505e-16;

static const Real unit_length = CONST_pc*1e3; // 1 kpc
static const Real unit_time   = CONST_yr*1e6; // 1 Myr
static const Real unit_density = CONST_amu;   // 1 mp/cm-3

static const Real unit_velocity = unit_length/unit_time; // in kpc/Myr

static const Real unit_q = (unit_density * pow(unit_velocity,3.0))/unit_length;

static const Real KELVIN = unit_velocity*unit_velocity*CONST_amu/CONST_kB;
*/


// Variable read from the input file

static int n1 = 64;
static int n2 = 64;
static int n3 = 64;

static int cooling_flag = 0;

static Real amb_rho = 1.0;

static Real cloud_radius;
static Real cloud_time;
static Real cloud_chi;

// Relevant temperatures
static Real T_floor   = 1e4;   // in K
static Real T_ceil    = 1e8;   // in K
static Real T_hot     = 1e7;
static Real T_hot_req = 1e7;
static Real T_cold    = 2*1e4;
static Real T_cut_mul = 0.6;
static Real T_cut     = T_cut_mul*T_hot_req;

// in terms of solar abundances
static Real Zsol = 1.0;
static Real Xsol = 1.0;

static Real X = Xsol * 0.7381;
static Real Z = Zsol * 0.0134;
static Real Y = 1 - X - Z;

static Real mu  = 1.0/(2.*X+ 3.*(1.-X-Z)/4.+ Z/2.);
static Real mue = 2.0/(1.0+X);
static Real muH = 1.0/X;

// Cloud position
static Real cloud_pos_x = 0.0;
static Real cloud_pos_y = 0.0;
static Real cloud_pos_z = 0.0;

// Magnetic fields
static Real B_x = 0.0;
static Real B_y = 0.0;
static Real B_z = 1.0;

// static Real mH = 1.0; //! Not required?

bool cooling_flag_print_count = false;

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

  T_floor      = pin->GetReal("problem","T_floor");
  T_ceil       = pin->GetReal("problem","T_ceil");
  T_hot        = pin->GetReal("problem","T_hot");
  T_hot_req    = pin->GetReal("problem","T_hot_req");
  T_cold       = pin->GetReal("problem","T_cold");
  T_cut_mul    = pin->GetReal("problem","T_cut_mul");
  T_cut        = pin->GetReal("problem","T_cut");

  // T_cut = T_cut_mul*T_hot_req;

  Xsol            = pin->GetReal("problem","Xsol");
  Zsol            = pin->GetReal("problem","Zsol");

  X = Xsol * 0.7381;
  Z = Zsol * 0.0134;
  Y = 1 - X - Z;

  mu  = 1.0/(2.*X+ 3.*(1.-X-Z)/4.+ Z/2.);
  mue = 2.0/(1.0+X);
  muH = 1.0/X;

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


void Cooling(MeshBlock *pmb, const Real time, const Real dt,
             const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
             AthenaArray<Real> &cons_scalar) {

  Real g = pmb->peos->GetGamma();

  Real sim_time = pmb->pmy_mesh->time;

  Real t_cloud = cloud_time;


  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    for (int j = pmb->js; j <= pmb->je; ++j) {
      for (int i = pmb->is; i <= pmb->ie; ++i) {

        Real temp = (prim(IPR,k,j,i) / prim(IDN,k,j,i)) * KELVIN * mu ;

        if (temp > T_floor) {

          Real *lam_para = Lam_file_read(temp);

          //printf("%f %f %f\n",lam_para[0],lam_para[1],lam_para[2]);


          Real temp_new = T_new(temp, lam_para[0], lam_para[1], lam_para[2],
                                prim(IDN,k,j,i), dt, 
                                mu, mue, muH, 
                                g, T_floor, T_ceil, T_cut);

                        // Defined in ../utils/townsend_cooling.hpp

          // ** NOTE: Above, dt is in code units to avoid overflow. unit_time is cancelled in 
          // ** calculation of T_new as we calculate (dt/tcool)

          cons(IEN,k,j,i) += ((temp_new-temp)/(KELVIN*mu))*prim(IDN,k,j,i)/(g-1);

          // ________________________________
          // FOR DEBUG PURPOSES
          if ((temp>T_cut) && ((temp_new-temp)!=0.0)){

            printf("T, delT: %lf %.60lf\n", temp, temp_new-temp);
          }
          // ________________________________


          
        }

        else{  // If T<=T_floor

          printf("+");
          cons(IEN,k,j,i) += ((T_floor-temp)/(KELVIN*mu))*prim(IDN,k,j,i)/(g-1);

        }

      }
    }
  }

  return;
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
  
    Cooling(pmb, time, dt,
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
  
  //* History outputs
  if (MAGNETIC_FIELDS_ENABLED) {

    AllocateUserHistoryOutput(8);

    EnrollUserHistoryOutput(0, rho_sum, "rho_sum");
    EnrollUserHistoryOutput(1, rho_sq_sum, "rho_sq_sum");
    EnrollUserHistoryOutput(2, c_s_sum, "c_s_sum");
    EnrollUserHistoryOutput(3, Pth_sum, "Pth_sum");
    EnrollUserHistoryOutput(4, PB_sum, "PB_sum");
    EnrollUserHistoryOutput(5, Bx_sum, "Bx_sum");
    EnrollUserHistoryOutput(6, By_sum, "By_sum");
    EnrollUserHistoryOutput(7, Bz_sum, "Bz_sum");
  }
  else {
    AllocateUserHistoryOutput(3);

    EnrollUserHistoryOutput(0, rho_sum, "rho_sum");
    EnrollUserHistoryOutput(1, rho_sq_sum, "rho_sq_sum");
    EnrollUserHistoryOutput(2, c_s_sum, "c_s_sum");
  }

  
  return;
}


//* To add clouds during restarts
void MeshBlock::InitUserMeshBlockData(ParameterInput *pin){

    Real g = peos->GetGamma();

    Real cloud_r  = cloud_radius;
    Real sim_time = pmy_mesh->time;
    Real t_cloud  = cloud_time;

    printf("InitUserMeshBlockData test: %lf \n", phydro->u(IDN,0,0,0));

    //* Rescale temperatures to T_hot_req
    if (temp_rescale_flag && (sim_time > 0.75*t_cloud)){

      Real local_T_sum_1 = 0.0;
      Real global_T_sum_1;

      // Calculate average temperature before rescaling
      for (int k = ks; k <= ke; ++k) {
        for (int j = js; j <= je; ++j) {
          for (int i = is; i <= ie; ++i) {

            Real temp  = (phydro->w(IPR,k,j,i) / phydro->u(IDN,k,j,i)) * KELVIN * mu ;
            
            local_T_sum_1 += temp;

          }
        }
      }

      MPI_Allreduce(&local_T_sum_1, &global_T_sum_1, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);

      Real N_cells = n1*n2*n3;
      Real T_avg = global_T_sum_1/N_cells; //* Average temperature

      // Rescaling the temperature
      for (int k = ks; k <= ke; ++k) {
        for (int j = js; j <= je; ++j) {
          for (int i = is; i <= ie; ++i) {

            Real temp  = (phydro->w(IPR,k,j,i) / phydro->u(IDN,k,j,i)) * KELVIN * mu ;
            Real T_new = temp*(T_hot_req/T_avg);

            phydro->u(IEN,k,j,i) += ((T_new-temp)/(KELVIN*mu))*phydro->u(IDN,k,j,i)/(g-1);

          }
        }
      }

     temp_rescale_flag = false; 

    }


    //* Inserting the cloud
    if (cloud_init_flag && (sim_time > 0.75*t_cloud)){

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
                    if(dist <= cloud_r*cloud_r){

                        Real temp = (phydro->w(IPR,k,j,i)/phydro->u(IDN,k,j,i)) * KELVIN * mu;

                        Real chi_lim = std::min(cloud_chi,temp/T_floor);

                        Real rho_temp = phydro->u(IDN,k,j,i);

                        //Initial energy ______________________________________//

                        Real TE_in = phydro->u(IEN,k,j,i);

                        Real KE_in = phydro->u(IM1,k,j,i)*phydro->u(IM1,k,j,i);
                        KE_in += phydro->u(IM2,k,j,i)*phydro->u(IM2,k,j,i);
                        KE_in += phydro->u(IM3,k,j,i)*phydro->u(IM3,k,j,i);
                        KE_in /= 2.0*phydro->u(IDN,k,j,i);

                        Real BE_in = 0.0;
                        if (MAGNETIC_FIELDS_ENABLED) {
                            BE_in += 0.5 * pfield->b.x1f(k,j,i) * pfield->b.x1f(k,j,i);
                            BE_in += 0.5 * pfield->b.x2f(k,j,i) * pfield->b.x2f(k,j,i);
                            BE_in += 0.5 * pfield->b.x3f(k,j,i) * pfield->b.x3f(k,j,i);
                        }

                        Real IE_in = TE_in - KE_in - BE_in;

                        //_______________________________________________________//

                        phydro->u(IDN,k,j,i) *= chi_lim;//*2.0;

                        phydro->u(IM1,k,j,i) /= sqrt(chi_lim);//*2.0;
                        phydro->u(IM2,k,j,i) /= sqrt(chi_lim);//*2.0;
                        phydro->u(IM3,k,j,i) /= sqrt(chi_lim);//*2.0;


                        //Final energy _________________________________________//

                        Real TE_fn = phydro->u(IEN,k,j,i);

                        Real KE_fn = phydro->u(IM1,k,j,i)*phydro->u(IM1,k,j,i);
                        KE_fn += phydro->u(IM2,k,j,i)*phydro->u(IM2,k,j,i);
                        KE_fn += phydro->u(IM3,k,j,i)*phydro->u(IM3,k,j,i);
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

                    }
                }
            }
        }

        printf("________________________________\n");
        printf("Cloud added!_______ %d, %lf /n",cloud_init_flag, sim_time);
        printf("________________________________\n");
        
        cloud_init_flag = false;

    }

    return;

}



//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  Real g   = pin->GetReal("hydro","gamma");


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


//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//  \brief
//========================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
}
