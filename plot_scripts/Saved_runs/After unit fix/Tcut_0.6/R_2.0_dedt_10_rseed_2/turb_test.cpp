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



// Variable read from the input file

static Real cloud_radius;
static Real cloud_time;
static Real cloud_chi;

static Real T_floor  = 1e4;  // in K
static Real T_ceil   = 1e8;   // in K
static Real T_hot = 1e7;
static Real T_cold = 2*1e4;
static Real T_cut_mul = 0.25;
static Real T_cut = T_cut_mul*T_hot;

static Real X = 1.0;
static Real Y = 0.0;
static Real Z = 0.0;

static Real cloud_pos_x = 0.0;
static Real cloud_pos_y = 0.0;
static Real cloud_pos_z = 0.0;

static Real mu  = 1.0/(2.*X+ 3.*(1.-X-Z)/4.+ Z/2.);
static Real mue = 2.0/(1.0+X);
static Real muH = 1.0/X;
 
static Real mH = 1.0;


Real* Lam_file_read (Real T){

  static Real *T_fit, *Lam_fit, *alp_fit;
  static int ntab;
  FILE *lam_read;

  if (T_fit == NULL){
    
    // printf("Reading power law fit file_____________\n");

    
    lam_read = fopen("../cooling_power_law/power_law_fit_Z_1.0.txt", "r");

    int line_num = 1;
    char ch;

    if(lam_read==NULL){
      perror("Error\n");
    }

    while (ch != EOF){
      ch = fgetc(lam_read);
      // printf("%c",ch);

      if (ch=='\n'){
        line_num = line_num+1;
      }
    }

    rewind(lam_read);

    printf("\nNumber of line: %d\n",line_num);


    // TRY TO FIND SOME OTHER WAY
    // CANNOT DELETE THE ALLOCATED MEMORY, AS IT IS NEEDED LATER!
    T_fit = new Real[line_num];
    Lam_fit = new Real[line_num];
    alp_fit = new Real[line_num];

    ntab = 0;
    while (fscanf(lam_read, "%lf %lf %lf\n", T_fit + ntab, Lam_fit + ntab, alp_fit+ntab)!=EOF) {
      // printf("Inside fscanf: %lf %lf %lf\n",*(T_fit+ntab),*(Lam_fit+ntab),*(alp_fit+ntab));
      ntab++;
    }

    fclose(lam_read);


  }

  int i_a = 0;
  int i_b = ntab-1;
  int mid;

  while (i_a!=i_b-1){

    mid = (i_a+i_b)/2;

    if (T<T_fit[mid]){
      i_b = mid;
    }
    else{
      i_a = mid;
    }

  }

  static Real return_arr[3];

  return_arr[0] = T_fit[i_a];
  return_arr[1] = Lam_fit[i_a];
  return_arr[2] = alp_fit[i_a];

  return return_arr;
}


Real Lambda(Real T, Real T_0, Real Lam0, Real alp, Real T_floor, Real T_ceil){
/*
Cooling function: returns electron cooling efficiency

Arguments:: Real T : Temperature
            Real T_floor : Floor temperature

Output:: Real Lam : Cooling funtion value, 0.0 if T<=T_floor  
*/

  Real Lam=0.0;

  if ((T >= T_floor) && (T < T_ceil)){     // >= must be used. Just > DOESN'T work
    // Lam = Lam0 + alp * std::log10(T/T_0);
    // Lam = pow(10,Lam0) * pow(T/T_0,alp);
    Lam = Lam0 * pow(T/T_0,alp);
  }

  // printf("Lam: %.50f\n",Lam);

  return Lam;
}

Real Y_TEF (Real T, Real T_0, Real Lam0, Real alp, Real T_floor){
/*
Dimensionless temporal evolution function (TEF)

Arguments:: Real T : Temperature
            Real T_floor : Floor temperature

Output:: Real Y_out
*/

  Real Y_out;

  if (alp==1.0) {
    Y_out = log(T_0/T);
  }
  else{
    Y_out = (1.0/(1.0-alp))*( 1.0- pow( T_0/T, alp-1.0 ) );
  }

  return Y_out;
}

Real Y_inv (Real y, Real T_0, Real Lam0, Real alp, Real T_floor){
/*
Inverse of Y_TEF function

Arguments:: Real y
            Real T_floor : Floor temperature

Output:: Real Y_inv_out: returns temperature corresponding y
*/

  Real Y_inv_out;

  if (alp==1.0){
    Y_inv_out = T_0*exp(-1.0*y);
  } 
  else{
    Y_inv_out = T_0*pow( 1.0-(1.0-alp)*y ,1.0/(1.0-alp));
  }

  return Y_inv_out;
}

Real tcool_calc(Real T_n, Real T_0, Real Lam0, Real alp, 
           Real rho, Real dt, Real mu, Real mue, Real muH, 
           Real gamma, Real T_floor, Real T_ceil){

  // Real k = CONST_kB;

  // printf("Lam: %.50f\n",CONST_kB);

  Real n_H = rho*unit_density/(muH*CONST_amu);

  Real lam_arr = Lambda(T_n,T_0,Lam0,alp,T_floor,T_ceil);

  Real p = rho*T_n/(KELVIN*mu);
  
  Real q = n_H*n_H*lam_arr/unit_q;

  Real tc = p/(q*(gamma - 1));

  return tc;  // in code units

  // Real Lam = Lambda(T_n,T_0,Lam0,alp,T_floor,T_ceil);

  // Real tcool = (CONST_kB*mue*muH*T_n)/((gamma-1.)*rho*mu*Lam);

  // return tcool/unit_time;
  // NOTE: tcool is returned in code units to avoid overflow
  // unit_time is cancelled in (dt/tcool)

}

Real T_new (Real T_n, Real T_0, Real Lam0, Real alp, 
            Real rho, Real dt, Real mu, Real mue, Real muH, 
            Real gamma, Real T_floor, Real T_ceil){
/*
Returns new temperture for each timestep

Arguments:: Real T_n : Initial temperature
            Real rho : Density
            Real dt  : Timestep
            Real mu  : Number density of particles
            Real mue : Electron number density
            Real muH : Proton number density
            Real gamma : Adiabatic index
            Real T_floor : Floor temperature

Output:: Real : returns final temperature after cooling
*/

  if (T_n >=T_cut){

    // printf("T_n: %.50f\n",T_n);

    return T_n;
  }

  Real T_nplus;

  // Real k = CONST_kB;

  Real tcool = tcool_calc(T_n, T_0, Lam0, alp, rho, dt, mu, mue, muH, gamma, T_floor, T_ceil);

  Real del_Y = (Lambda(T_0,T_0,Lam0,alp,T_floor,T_ceil)/Lambda(T_n,T_0,Lam0,alp,T_floor,T_ceil));

  del_Y *= (T_n/T_0) * (dt/tcool);

  // printf("T_cut: %.50f\n",T_cut);

  T_nplus = Y_inv( Y_TEF(T_n,T_0,Lam0,alp,T_floor) + del_Y,T_0, Lam0, alp, T_floor );


  if (std::isnan(T_nplus)){
    return T_floor;
  }
  else if (T_nplus>T_ceil){
    return T_ceil;
  } 
  else if (T_nplus<T_floor){
    return T_floor;
  }
  else{
    return T_nplus;
  }


}


void Cooling(MeshBlock *pmb, const Real time, const Real dt,
             const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
             AthenaArray<Real> &cons_scalar) {

  Real g = pmb->peos->GetGamma();

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
                                g, T_floor, T_ceil);

          // NOTE: Above, dt is in code units to avoid overflow. unit_time is cancelled in 
          // calculation of T_new as we calculate (dt/tcool)

          //if (std::abs(temp_new-temp)>=0.1){
          //  printf("del_T: %f\n",temp_new-temp);
          //}

          if ((temp>T_cut) && ((temp_new-temp)!=0.0)):

            printf("T, delT: %lf %lf\n", temp, temp_new-temp);

          cons(IEN,k,j,i) += ((temp_new-temp)/(KELVIN*mu))*prim(IDN,k,j,i)/(g-1);
        }
      }
    }
  }

  return;
}


void Cloud_init(MeshBlock *pmb, const Real time, const Real dt,
             const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
             AthenaArray<Real> &cons_scalar) {

  Real cloud_r = cloud_radius;
  
  Real sim_time = pmb->pmy_mesh->time;
  Real sim_cycle = pmb->pmy_mesh->ncycle;

  Real g = pmb->peos->GetGamma();

  static int count;  
  static int ncyle;

  Real t_cloud = cloud_time;

  if (count==0){

    if ((sim_time<t_cloud+dt) && (sim_time>t_cloud-dt)){
      
      for (int k = pmb->ks; k <= pmb->ke; ++k) {
      Real z_coord = pmb->pcoord->x3v(k);

        for (int j = pmb->js; j <= pmb->je; ++j) {
          Real y_coord = pmb->pcoord->x2v(j);

          for (int i = pmb->is; i <= pmb->ie; ++i) {
            Real x_coord = pmb->pcoord->x1v(i);
            
            Real dist = (x_coord-cloud_pos_x)*(x_coord-cloud_pos_x);
            dist += (y_coord-cloud_pos_y)*(y_coord-cloud_pos_y);
            dist += (z_coord-cloud_pos_z)*(z_coord-cloud_pos_z);

            if(dist <= cloud_r*cloud_r){

              Real temp = (prim(IPR,k,j,i)/cons(IDN,k,j,i)) * KELVIN * mu;

              Real chi_lim = std::min(cloud_chi,temp/T_floor);

              cons(IDN,k,j,i) *= chi_lim;

            }
            
          }
        }
      }

      printf("________________________________\n");
      printf("Cloud added!\n");
      printf("________________________________\n");
      
      count += 1;
    }
  
  }

  // printf("%d \n", count);

  
  return;
}


void Source(MeshBlock *pmb, const Real time, const Real dt,
             const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
             AthenaArray<Real> &cons_scalar) {

	Cloud_init(pmb, time, dt,
             prim, prim_scalar,
             bcc, cons,
             cons_scalar);

  Cooling(pmb, time, dt,
             prim, prim_scalar,
             bcc, cons,
             cons_scalar);

  return;

}


Real cold_gas(MeshBlock *pmb, int iout)
{
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


Real rho_sum(MeshBlock *pmb, int iout)
{
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


Real rho_sq_sum(MeshBlock *pmb, int iout)
{
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

  cloud_radius = pin->GetReal("problem","cloud_radius");
  cloud_time   = pin->GetReal("problem","cloud_time");
  cloud_chi    = pin->GetReal("problem","cloud_chi");

  T_floor      = pin->GetReal("problem","T_floor");
  T_ceil       = pin->GetReal("problem","T_ceil");
  T_hot        = pin->GetReal("problem","T_hot");
  T_cold        = pin->GetReal("problem","T_cold");
  T_cut_mul    = pin->GetReal("problem","T_cut_mul");

  T_cut = T_cut_mul*T_hot;

  X            = pin->GetReal("problem","X");
  Y            = pin->GetReal("problem","Y");
  Z            = pin->GetReal("problem","Z");

  cloud_pos_x  = pin->GetReal("problem","cloud_pos_x");
  cloud_pos_y  = pin->GetReal("problem","cloud_pos_y");
  cloud_pos_z  = pin->GetReal("problem","cloud_pos_z");


  //====================  NEW  ======================================

  EnrollUserExplicitSourceFunction(Source);
  printf("\n__________Source! Inside InitUserMeshData________________\n");

  AllocateUserHistoryOutput(3);

  EnrollUserHistoryOutput(0, cold_gas, "cold_gas");
  EnrollUserHistoryOutput(1, rho_sum, "rho_sum");
  EnrollUserHistoryOutput(2, rho_sq_sum, "rho_sq_sum");


  //=================================================================

  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  Real g = pin->GetReal("hydro","gamma");

  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        phydro->u(IDN,k,j,i) = 1.0;

        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;

        if (NON_BAROTROPIC_EOS) {
          phydro->u(IEN,k,j,i) = (T_hot/(KELVIN*mu))*phydro->u(IDN,k,j,i)/(g-1);
        }
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
