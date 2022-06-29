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

#include <algorithm>  // max()
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

#include "../inputs/hdf5_reader.hpp"  // HDF5ReadRealArray()

#include "../utils/townsend_cooling.hpp"

#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

// User headers
#include <fstream>    // for file io
#include <iomanip>    // for std::setprecision


// * Defined in townsend_cooling.hpp
//TODO: Maybe move the units to another file
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

static Real cloud_pos_x = 0.0;
static Real cloud_pos_y = 0.0;
static Real cloud_pos_z = 0.0;


static Real B_x = 0.0;
static Real B_y = 0.0;
static Real B_z = 1.0;



// static Real mH = 1.0; //! Not required?

static Real T_avg_before_rescale = 0.0;


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

  Xsol            = pin->GetReal("problem","X");
  Zsol            = pin->GetReal("problem","Z");

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

          // printf("!! Temperature floor applied ___________ !!");
          // printf("T, T_floor, prs: %lf %lf %lf \n", temp, T_floor, prim(IPR,k,j,i));
          // printf("!! _____________________________________ !!");
          printf("+");
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



//TODO: Move history output functions to another file

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
    if (pmb->pmy_mesh->time < cloud_time){

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


  // Read static variables from the input file
  read_input(pin);


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
  // Determine locations of initial values
  std::string input_filename = pin->GetString("problem", "input_filename");
  std::string dataset_cons = pin->GetString("problem", "dataset_cons");
  int index_dens = pin->GetInteger("problem", "index_dens");
  int index_mom1 = pin->GetInteger("problem", "index_mom1");
  int index_mom2 = pin->GetInteger("problem", "index_mom2");
  int index_mom3 = pin->GetInteger("problem", "index_mom3");
  int index_etot = pin->GetInteger("problem", "index_etot");
  std::string dataset_b1 = pin->GetString("problem", "dataset_b1");
  std::string dataset_b2 = pin->GetString("problem", "dataset_b2");
  std::string dataset_b3 = pin->GetString("problem", "dataset_b3");

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

  int start_field_mem[3];
  start_field_mem[0] = ks;
  start_field_mem[1] = js;
  start_field_mem[2] = is;

  int count_field_mem[3];

  // Set magnetic field values from file
  if (MAGNETIC_FIELDS_ENABLED) {
    // Set B1
    count_field_file[1] = block_size.nx3;
    count_field_file[2] = block_size.nx2;
    count_field_file[3] = block_size.nx1 + 1;
    count_field_mem[0] = block_size.nx3;
    count_field_mem[1] = block_size.nx2;
    count_field_mem[2] = block_size.nx1 + 1;
    HDF5ReadRealArray(input_filename.c_str(), dataset_b1.c_str(), 4, start_field_file,
                      count_field_file, 3, start_field_mem,
                      count_field_mem, pfield->b.x1f, true);

    // Set B2
    count_field_file[1] = block_size.nx3;
    count_field_file[2] = block_size.nx2 + 1;
    count_field_file[3] = block_size.nx1;
    count_field_mem[0] = block_size.nx3;
    count_field_mem[1] = block_size.nx2 + 1;
    count_field_mem[2] = block_size.nx1;
    HDF5ReadRealArray(input_filename.c_str(), dataset_b2.c_str(), 4, start_field_file,
                      count_field_file, 3, start_field_mem,
                      count_field_mem, pfield->b.x2f, true);

    // Set B3
    count_field_file[1] = block_size.nx3 + 1;
    count_field_file[2] = block_size.nx2;
    count_field_file[3] = block_size.nx1;
    count_field_mem[0] = block_size.nx3 + 1;
    count_field_mem[1] = block_size.nx2;
    count_field_mem[2] = block_size.nx1;
    HDF5ReadRealArray(input_filename.c_str(), dataset_b3.c_str(), 4, start_field_file,
                      count_field_file, 3, start_field_mem,
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
          HDF5ReadRealArray(input_filename.c_str(), dataset_b1.c_str(), 4,
                            start_field_file, count_field_file, 3,
                            start_field_mem, count_field_mem,
                            pfield->b.x1f, true, true);

          count_field_file[1] = block_size.nx3;
          count_field_file[2] = block_size.nx2 + 1;
          count_field_file[3] = block_size.nx1;
          count_field_mem[0] = block_size.nx3;
          count_field_mem[1] = block_size.nx2 + 1;
          count_field_mem[2] = block_size.nx1;
          HDF5ReadRealArray(input_filename.c_str(), dataset_b2.c_str(), 4,
                            start_field_file, count_field_file, 3,
                            start_field_mem, count_field_mem,
                            pfield->b.x2f, true, true);

          count_field_file[1] = block_size.nx3 + 1;
          count_field_file[2] = block_size.nx2;
          count_field_file[3] = block_size.nx1;
          count_field_mem[0] = block_size.nx3 + 1;
          count_field_mem[1] = block_size.nx2;
          count_field_mem[2] = block_size.nx1;
          HDF5ReadRealArray(input_filename.c_str(), dataset_b3.c_str(), 4,
                            start_field_file, count_field_file, 3,
                            start_field_mem, count_field_mem,
                            pfield->b.x3f, true, true);
                            
        }
      }
    }
  }
#endif
  return;
}



// * ! First do the the simple restart test
/*

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  Real g   = pin->GetReal("hydro","gamma");


  // Read static variables from the input file
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
*/

//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//  \brief
//========================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
}
