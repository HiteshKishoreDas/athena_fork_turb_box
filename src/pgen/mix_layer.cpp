//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mix_layer.cpp
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

// #include "../utils/townsend_cooling.hpp"  // T_new()
#include "../utils/townsend_cooling_max.hpp" // T_new()

#include "../utils/hst_func.hpp"             // All history output functions
#include "../utils/code_units.hpp"           // Code units and constants

//* ___________________________
//* For Max's townsend cooling

using namespace std;

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

static Real L1 = 0.0;
static Real L2 = 0.0;
static Real L3 = 0.0;

static int cooling_flag = 1.0;

static Real amb_rho = 1.0;

static Real front_thick = 2.5;
static Real v_shear = 100 * (1.023*1e-3); // 100 km/s in kpc/Myr

static Real knx_KH  = 1.0;
static Real kny_KH  = 1.0; 
static Real amp_KH = 0.01;    // Amplitude = amp_KH * v_shear

static Real B_x = 0.0;
static Real B_y = 0.0;
static Real B_z = 1.0;

static bool cooling_flag_print_count = false;
static bool DEBUG_FLAG = false;


void read_input (ParameterInput *pin){
  /*
  Reads all the static variables from the input file
  */

  L1 = pin->GetReal("mesh","x1max") - pin->GetReal("mesh","x1min");
  L2 = pin->GetReal("mesh","x2max") - pin->GetReal("mesh","x2min");
  L3 = pin->GetReal("mesh","x3max") - pin->GetReal("mesh","x3min");

  cooling_flag    = pin->GetInteger("problem","cooling_flag");

  amb_rho      = pin->GetReal("problem","amb_rho");

  front_thick  = pin->GetReal("problem","front_thickness");
  v_shear      = pin->GetReal("problem","v_shear");

  knx_KH = pin->GetReal("problem","knx_KH");
  kny_KH = pin->GetReal("problem","kny_KH");
  amp_KH = pin->GetReal("problem","amp_KH");

  T_floor      = pin->GetReal("problem","T_floor");
  T_ceil       = pin->GetReal("problem","T_ceil");
  T_hot        = pin->GetReal("problem","T_hot");
  T_cold       = pin->GetReal("problem","T_cold");
  T_cut_mul    = pin->GetReal("problem","T_cut_mul");

  T_cut = T_cut_mul*T_hot;

  Xsol            = pin->GetReal("problem","Xsol");
  Zsol            = pin->GetReal("problem","Zsol");

  X = Xsol * 0.7381;
  Z = Zsol * 0.0134;
  Y = 1 - X - Z;

  mu  = 1.0/(2.*X+ 3.*(1.-X-Z)/4.+ Z/2.);
  mue = 2.0/(1.0+X);
  muH = 1.0/X;

  Lambda_fac = pin->GetReal("problem","Lambda_fac");
  DEBUG_FLAG = pin->GetReal("problem","DEBUG_FLAG");

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

  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    for (int j = pmb->js; j <= pmb->je; ++j) {
      for (int i = pmb->is; i <= pmb->ie; ++i) {

        Real temp = (prim(IPR,k,j,i) / cons(IDN,k,j,i)) * KELVIN * mu ;

        //TODO: Move NaN checks to Cooling header file
        if (std::isnan(temp)) {

          printf("[cooling] nan detected for T at i,j,k : (%d, %d, %d) \n", i,j,k); 
          printf("temp = %lf \n", temp ); 
          printf("rho  = %lf \n", cons(IDN,k,j,i)); 
          printf("IM1  = %lf \n", cons(IM1,k,j,i)); 
          printf("IM2  = %lf \n", cons(IM2,k,j,i)); 
          printf("IM3  = %lf \n", cons(IM3,k,j,i)); 

          Real E_kin = cons(IM1,k,j,i)*cons(IM1,k,j,i);
          E_kin     += cons(IM2,k,j,i)*cons(IM2,k,j,i);
          E_kin     += cons(IM3,k,j,i)*cons(IM3,k,j,i);
          E_kin     /= 2*cons(IDN,k,j,i);

          Real E_mag = 0.0;

          if (MAGNETIC_FIELDS_ENABLED) {
            E_mag += prim(IB1,k,j,i)*prim(IB1,k,j,i);
            E_mag += prim(IB2,k,j,i)*prim(IB2,k,j,i);
            E_mag += prim(IB3,k,j,i)*prim(IB3,k,j,i);
            E_mag *= 0.5;
          }

          // Set current temperature to T_floor
          cons(IEN,k,j,i) = E_kin + E_mag + (T_floor/(KELVIN*mu))*cons(IDN,k,j,i)/(g-1);
          continue; 

        } // End of NaN check on current temperature

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
            Real cLfac    = Lambda_fac;

            Real temp_new = max(cooler->townsend(temp_cgs,rho_cgs,dt_cgs, cLfac), T_floor);

            if  (std::isnan(temp_new)) {

              printf("[cooling] nan detected for T_new at i,j,k : (%d, %d, %d) \n", i,j,k); 
              printf("temp = %lf \n", temp_cgs ); 
              printf("rho  = %lf \n", rho); 
              printf("IM1  = %lf \n", cons(IM1,k,j,i)); 
              printf("IM2  = %lf \n", cons(IM2,k,j,i)); 
              printf("IM3  = %lf \n", cons(IM3,k,j,i)); 

              // temp_new set to T_floor
              Real ccool = ((T_floor-temp)/(KELVIN*mu))*cons(IDN,k,j,i)/(g-1); 
              cons(IEN,k,j,i) += ccool;

            } //* End of NaN check
            
            else {
            
              Real ccool = ((temp_new-temp)/(KELVIN*mu))*cons(IDN,k,j,i)/(g-1);

              //TODO: Check the volume check with homogenous box
              //TODO: Put a test with ccool = 1

              cons(IEN,k,j,i) += ccool;
              
              if (DEBUG_FLAG){
                total_cooling -= 1.0;
              }
              else{
                total_cooling -= ccool;
              }


            } // End of else for NaN check

          } // End of T<T_cut check
        } // End of T>T_floor check

        else{  // If T<=T_floor

          printf("+");

          Real ccool = ((T_floor-temp)/(KELVIN*mu))*cons(IDN,k,j,i)/(g-1); 
          cons(IEN,k,j,i) += ccool;
          // total_cooling -= ccool;

        } // End of else, T<=T_floor

      }
    }
  } // End of loop over Meshblock

  return;
}


Real hst_total_cooling(MeshBlock *pmb, int iout) {
  if(pmb->lid == 0)
    return total_cooling;
  else
    return 0;
}


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


  //* Read all the static variables from the input file
  read_input(pin);

  //* Enroll the Source terms
  EnrollUserExplicitSourceFunction(Source);

  cooler = make_unique<Cooling>();

  //* History outputs
  if (MAGNETIC_FIELDS_ENABLED) {

    AllocateUserHistoryOutput(10);

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
  }
  else {

    AllocateUserHistoryOutput(5);

    EnrollUserHistoryOutput(0, rho_sum, "rho_sum");
    EnrollUserHistoryOutput(1, rho_sq_sum, "rho_sq_sum");
    EnrollUserHistoryOutput(2, c_s_sum, "c_s_sum");
    EnrollUserHistoryOutput(3, cold_gas, "cold_gas");
    EnrollUserHistoryOutput(4, hst_total_cooling, "total_cooling");

  }

  //=================================================================

  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  Real g   = pin->GetReal("hydro","gamma");

  //* Read all the static variables from the input file
  read_input(pin);

  Real A_KH = amp_KH * v_shear;
  Real k_x  = 2*PI*knx_KH/L1;
  Real k_y  = 2*PI*kny_KH/L2;

  printf("k_x, k_y: %lf \t %lf \n", k_x, k_y);
  printf("L1, L2, L3: %lf \t %lf \t %lf \n", L1, L2, L3);

  Real rho_cold = amb_rho * T_hot/T_floor;

  Real x1, x2, x3;
  Real x1f, x2f, x3f;

  for (int k=ks; k<=ke; k++) {

    x3 = pcoord->x3v(k);
    x3f = pcoord->x3f(k);

    for (int j=js; j<=je; j++) {

      x2 = pcoord->x2v(j);
      x2f = pcoord->x2f(j);

      for (int i=is; i<=ie; i++) {

        x1 = pcoord->x1v(i);
        x1f = pcoord->x1f(i);   
        
        phydro->u(IDN,k,j,i) = rho_cold + (0.5*(amb_rho - rho_cold)) * (1 + tanh(x3/front_thick));

        phydro->u(IM1,k,j,i) = phydro->u(IDN,k,j,i) * (v_shear/2) * (1 + tanh(x3/front_thick));
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = phydro->u(IDN,k,j,i) * A_KH * exp(-1.0*x3*x3/front_thick/front_thick) * sin(k_x*x1) * sin(k_y*x2);

        if (NON_BAROTROPIC_EOS) {
          // Kinetic energy
          phydro->u(IEN,k,j,i)  = phydro->u(IM1,k,j,i)*phydro->u(IM1,k,j,i);
          phydro->u(IEN,k,j,i) += phydro->u(IM2,k,j,i)*phydro->u(IM2,k,j,i);
          phydro->u(IEN,k,j,i) += phydro->u(IM3,k,j,i)*phydro->u(IM3,k,j,i);
          phydro->u(IEN,k,j,i) /= 2.0*phydro->u(IDN,k,j,i);

          //Internal energy
          phydro->u(IEN,k,j,i) += (T_hot/(KELVIN*mu))*amb_rho/(g-1);
        }

        // printf("rho: %f\n",phydro->u(IDN,k,j,i));
      }
    }
  }

  // initialize interface B
  if (MAGNETIC_FIELDS_ENABLED) {
    
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

            // Magnetic energy
            phydro->u(IEN,k,j,i) += 0.5*(SQR(B_x) + SQR(B_y) + SQR(B_z));
          }
        }
      }
    }
  }

  printf("________________________________________________\n");
  printf("Magnetic fields initialised!  \n");
  printf("________________________________________________\n");

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

        Real lum_cell = 0.0;
        Real temp = (phydro->u(IPR,k,j,i) / phydro->u(IDN,k,j,i)) * KELVIN * mu ;

        if (temp > T_floor) {

            //* For Max's townsend cooling

            if (temp<T_cut){
              Real rho      = phydro->u(IDN,k,j,i);
              Real temp_cgs = temp;
              Real rho_cgs  = rho * unit_density;
              Real dt_cgs   = pmy_mesh->dt * unit_time;
              Real cLfac    = 1.0;

              Real temp_new = max(cooler->townsend(temp_cgs,rho_cgs,dt_cgs, 1.0), T_floor);

              Real ccool_2 = ((temp_new-temp)/(KELVIN*mu))*phydro->u(IDN,k,j,i)/(g-1);

              lum_cell-= ccool_2;
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
