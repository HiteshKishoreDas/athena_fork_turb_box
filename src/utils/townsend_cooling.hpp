#ifndef TOWNSEND_H
#define TOWNSEND_H

// C++ headers
#include <cmath>
#include <ctime>
#include <sstream>
#include <stdexcept>

// User headers
#include <fstream>    // for file io
#include <iomanip>    // for std::setprecision

// Athena headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"


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



Real* Lam_file_read (Real T){

  static Real *T_fit, *Lam_fit, *alp_fit;
  static int ntab;
  FILE *lam_read;

  if (T_fit == NULL){
    
    printf("Reading power law fit file_____________\n");

    lam_read = fopen("../../cooling_power_law/power_law_fit_Z_1.0.txt", "r");

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
            Real gamma, Real T_floor, Real T_ceil,
            Real T_cut){
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


#endif