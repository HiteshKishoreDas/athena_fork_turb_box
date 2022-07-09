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

#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

// User headers
#include <fstream>    // for file io
#include <iomanip>    // for std::setprecision

#include "../utils/townsend_cooling.hpp"

// Variable read from the input file

static Real L1 = 0.0;
static Real L2 = 0.0;
static Real L3 = 0.0;

static Real gamma_in = 5.0/3.0;

static int cooling_flag = 0;

static Real amb_rho = 1.0;

static Real front_thick = 2.5;
static Real v_shear = 100 * (1.023*1e-3); // 100 km/s in kpc/Myr

static Real knx_KH  = 1.0;
static Real kny_KH  = 1.0; 
static Real amp_KH = 0.01;    // Amplitude = amp_KH * v_shear

static Real T_floor  = 1e4;  // in K
static Real T_ceil   = 1e8;   // in K
static Real T_hot = 1e7;
static Real T_cold = 2*1e4;
static Real T_cut_mul = 0.6;
static Real T_cut;     // = T_cut_mul*T_hot;

static Real X = 1.0;
static Real Y = 0.0;
static Real Z = 0.0;

static Real B_x = 0.0;
static Real B_y = 0.0;
static Real B_z = 1.0;


static Real mu  = 1.0/(2.*X+ 3.*(1.-X-Z)/4.+ Z/2.);
static Real mue = 2.0/(1.0+X);
static Real muH = 1.0/X;
 
static Real mH = 1.0;



void read_input (ParameterInput *pin){
  /*
  Reads all the static variables from the input file
  */

  L1 = pin->GetReal("mesh","x1max") - pin->GetReal("mesh","x1min");
  L2 = pin->GetReal("mesh","x2max") - pin->GetReal("mesh","x2min");
  L3 = pin->GetReal("mesh","x3max") - pin->GetReal("mesh","x3min");

  gamma_in   = pin->GetReal("hydro","gamma");

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

  X            = pin->GetReal("problem","X");
  Y            = pin->GetReal("problem","Y");
  Z            = pin->GetReal("problem","Z");

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

  // printf("T_cut in Cooling: %lf \n", T_cut);

  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    for (int j = pmb->js; j <= pmb->je; ++j) {
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real temp = (prim(IPR,k,j,i) / prim(IDN,k,j,i)) * KELVIN * mu ;

        if (temp > T_floor) {

          if (cooling_flag){

            Real *lam_para = Lam_file_read(temp);

            //printf("%f %f %f\n",lam_para[0],lam_para[1],lam_para[2]);


            Real temp_new = T_new(temp, lam_para[0], lam_para[1], lam_para[2],
                                  prim(IDN,k,j,i), dt, 
                                  mu, mue, muH, 
                                  g, T_floor, T_ceil);

            // NOTE: Above, dt is in code units to avoid overflow. unit_time is cancelled in 
            // calculation of T_new as we calculate (dt/tcool)

            cons(IEN,k,j,i) += ((temp_new-temp)/(KELVIN*mu))*prim(IDN,k,j,i)/(g-1);

            // printf("T, delT: %lf %.60lf\n", temp, temp_new-temp);

            // ________________________________
            // FOR DEBUG PURPOSES
            if ((temp>T_cut) && ((temp_new-temp)!=0.0)){
              printf("T, delT: %lf %.60lf\n", temp, temp_new-temp);
            }
            // ________________________________
            
          }

        }

        else if (temp < T_floor){

          cons(IEN,k,j,i) += ((T_floor-temp)/(KELVIN*mu))*prim(IDN,k,j,i)/(g-1);

        }

      }
    }
  }

  return;
}


static int cooling_flag_print_count = 0;

void Source(MeshBlock *pmb, const Real time, const Real dt,
             const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
             AthenaArray<Real> &cons_scalar) {

  if (cooling_flag!=0){

    if (cooling_flag_print_count==0){ 
      printf("___________________________________________\n");
      printf("!! Cooling included .......................\n");
      printf("___________________________________________\n");

      cooling_flag_print_count += 1;
    }

    Cooling(pmb, time, dt,
            prim, prim_scalar,
            bcc, cons,
            cons_scalar);
  }

  return;

}


Real cold_gas(MeshBlock *pmb, int iout){

  Real cold_gas_mass=0;

  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  
  for(int k=ks; k<=ke; k++) {
    for(int j=js; j<=je; j++) {
      for(int i=is; i<=ie; i++) {
        
        Real rho = pmb->phydro->u(IDN,k,j,i);
        Real prs = pmb->phydro->w(IPR,k,j,i);

        Real temp = (prs / rho) * KELVIN * mu ;

        if (temp <= T_cold){
          cold_gas_mass += rho*pmb->pcoord->GetCellVolume(k,j,i);
        }

      }
    }
  }

  return cold_gas_mass;
}

Real rho_sum(MeshBlock *pmb, int iout){

  Real rho_sum    = 0;

  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  
  for(int k=ks; k<=ke; k++) {
    for(int j=js; j<=je; j++) {
      for(int i=is; i<=ie; i++) {
        
        Real rho = pmb->phydro->u(IDN,k,j,i);

        rho_sum += rho;

      }
    }
  }

  return rho_sum;
}

Real rho_sq_sum(MeshBlock *pmb, int iout){

  Real rho_sq_sum = 0;

  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  
  for(int k=ks; k<=ke; k++) {
    for(int j=js; j<=je; j++) {
      for(int i=is; i<=ie; i++) {
        
        Real rho = pmb->phydro->u(IDN,k,j,i);
        rho_sq_sum += rho*rho;

      }
    }
  }

  return rho_sq_sum;
}

Real c_s_sum(MeshBlock *pmb, int iout){

  Real c_s_sum = 0;
  Real gamma = pmb->peos->GetGamma();

  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;

  for(int k=ks; k<=ke; k++) {
    for(int j=js; j<=je; j++) {
      for(int i=is; i<=ie; i++) {
        
        Real rho = pmb->phydro->u(IDN,k,j,i);
        Real prs = pmb->phydro->w(IPR,k,j,i);

        c_s_sum += sqrt(gamma*prs/rho);

      }
    }
  }

  return c_s_sum;
}

Real tcool_sum(MeshBlock *pmb, int iout){

  Real tcool_sum = 0;
  
  if (cooling_flag!=0){
    

    Real gamma = pmb->peos->GetGamma();
    int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;

    Real dt = pmb->pmy_mesh->dt;

    for(int k=ks; k<=ke; k++) {
      for(int j=js; j<=je; j++) {
        for(int i=is; i<=ie; i++) {

          Real rho = pmb->phydro->u(IDN,k,j,i);
          Real prs = pmb->phydro->w(IPR,k,j,i);
          Real temp = (prs / rho) * KELVIN * mu ;

          Real *lam_para = Lam_file_read(temp);

          Real T_0  = lam_para[0];
          Real Lam0 = lam_para[1];
          Real alp  = lam_para[2];
      
          Real tcool = tcool_calc(temp, T_0, Lam0, alp, rho, dt, mu, mue, muH, gamma, T_floor, T_ceil);

          tcool_sum += tcool;

        }
      }
    }
  
  }

  return tcool_sum;

}

Real Pth_sum(MeshBlock *pmb, int iout){

  Real Pth_sum = 0;
  

  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;


  for(int k=ks; k<=ke; k++) {
    for(int j=js; j<=je; j++) {
      for(int i=is; i<=ie; i++) {

        Pth_sum += pmb->phydro->w(IPR,k,j,i);
        
      }
    }
  }
    
  

  return Pth_sum;

}

Real PB_sum(MeshBlock *pmb, int iout){

  Real PB_sum = 0;
  
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;  

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        PB_sum += pmb->pfield->b.x1f(k,j,i)*pmb->pfield->b.x1f(k,j,i);
        PB_sum += pmb->pfield->b.x2f(k,j,i)*pmb->pfield->b.x2f(k,j,i);
        PB_sum += pmb->pfield->b.x3f(k,j,i)*pmb->pfield->b.x3f(k,j,i);
      }
    }
  }

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      PB_sum += pmb->pfield->b.x1f(k,j,ie+1)*pmb->pfield->b.x1f(k,j,ie+1);
    }
  }

  for (int k=ks; k<=ke; ++k) {
    for (int i=is; i<=ie; ++i) {
      PB_sum += pmb->pfield->b.x2f(k,je+1,i)*pmb->pfield->b.x2f(k,je+1,i);
    }
  }

  for (int j=js; j<=je; ++j) {
    for (int i=is; i<=ie; ++i) {
      PB_sum += pmb->pfield->b.x3f(ke+1,j,i)*pmb->pfield->b.x3f(ke+1,j,i);
    }
  }

  return PB_sum/2.0;

}

Real Bx_sum(MeshBlock *pmb, int iout){

  Real Bx_sum = 0;
  
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;  

  for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie+1; ++i) {
          Bx_sum += pmb->pfield->b.x1f(k,j,i);
        }
      }
    }

  return Bx_sum;

}

Real By_sum(MeshBlock *pmb, int iout){

  Real By_sum = 0;
  
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;

  for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je+1; ++j) {
        for (int i=is; i<=ie; ++i) {
          By_sum += pmb->pfield->b.x2f(k,j,i);
        }
      }
    }

  return By_sum;

}

Real Bz_sum(MeshBlock *pmb, int iout){

  Real Bz_sum = 0;
  
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;

  for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          Bz_sum += pmb->pfield->b.x3f(k,j,i);
        }
      }
    }

  return Bz_sum;

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
//   turb_flag = pin->GetInteger("problem","turb_flag");
//   if (turb_flag != 0) {
// #ifndef FFT
//     std::stringstream msg;
//     msg << "### FATAL ERROR in TurbulenceDriver::TurbulenceDriver" << std::endl
//         << "non zero Turbulence flag is set without FFT!" << std::endl;
//     ATHENA_ERROR(msg);
//     return;
// #endif
//   }


  // Read all the static variables from the input file
  read_input(pin);

  // printf("T_cut in InitUserMeshData: %lf \n", T_cut);

  //====================  NEW  ======================================


  EnrollUserExplicitSourceFunction(Source);
  printf("\n__________Source! Inside InitUserMeshData________________\n");


  if (MAGNETIC_FIELDS_ENABLED) {

    AllocateUserHistoryOutput(10);

    EnrollUserHistoryOutput(0, cold_gas, "cold_gas");
    EnrollUserHistoryOutput(1, rho_sum, "rho_sum");
    EnrollUserHistoryOutput(2, rho_sq_sum, "rho_sq_sum");
    EnrollUserHistoryOutput(3, c_s_sum, "c_s_sum");
    EnrollUserHistoryOutput(4, tcool_sum, "tcool_sum");
    EnrollUserHistoryOutput(5, Pth_sum, "Pth_sum");
    EnrollUserHistoryOutput(6, PB_sum, "PB_sum");
    EnrollUserHistoryOutput(7, Bx_sum, "Bx_sum");
    EnrollUserHistoryOutput(8, By_sum, "By_sum");
    EnrollUserHistoryOutput(9, Bz_sum, "Bz_sum");

  }
  else {

    AllocateUserHistoryOutput(5);

    EnrollUserHistoryOutput(0, cold_gas, "cold_gas");
    EnrollUserHistoryOutput(1, rho_sum, "rho_sum");
    EnrollUserHistoryOutput(2, rho_sq_sum, "rho_sq_sum");
    EnrollUserHistoryOutput(3, c_s_sum, "c_s_sum");
    EnrollUserHistoryOutput(4, tcool_sum, "tcool_sum");

  }

  //=================================================================

  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  // Read all the static variables from the input file
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
          phydro->u(IEN,k,j,i) += (T_hot/(KELVIN*mu))*amb_rho/(gamma_in-1);
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

//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//  \brief
//========================================================================================

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin){
  AllocateUserOutputVariables(1);
  return;
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin){

  for(int k=ks; k<=ke; k++) {
    for(int j=js; j<=je; j++) {
      for(int i=is; i<=ie; i++) {

        Real temp = ( phydro->u(IPR,k,j,i) / phydro->u(IDN,k,j,i) ) * KELVIN * mu ;
        Real *lam_para = Lam_file_read(temp);

        Real temp_new = T_new(temp, lam_para[0], lam_para[1], lam_para[2],
                                  phydro->u(IDN,k,j,i), pmy_mesh->dt, 
                                  mu, mue, muH, 
                                  gamma_in, T_floor, T_ceil);


        // Real tcool = tcool_calc(temp, lam_para[0], lam_para[1], lam_para[2],
        //             phydro->u(IDN,k,j,i), 1.0, 
        //             mu, mue, muH, 
        //             gamma_in, T_floor, T_ceil);

        
        Real dE = ((temp_new-temp)/(KELVIN*mu))*phydro->u(IDN,k,j,i)/(gamma_in-1);

        user_out_var(0,k,j,i) = dE/pmy_mesh->dt;      
      
      }
    }
  }
}

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
}
