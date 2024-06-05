#ifndef CODE_UNITS_H
#define CODE_UNITS_H

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

// in terms of solar abundances
static Real Zsol = 1.0;
static Real Xsol = 1.0;

static Real X = Xsol * 0.7381;
static Real Z = Zsol * 0.0134;
static Real Y = 1 - X - Z;

static Real mu  = 1.0/(2.*X+ 3.*(1.-X-Z)/4.+ Z/2.);
static Real mue = 2.0/(1.0+X);
static Real muH = 1.0/X;

// Relevant temperatures
static Real T_floor;   // in K
static Real T_ceil;   // in K
static Real T_hot;
static Real T_hot_req;
static Real T_cold;
static Real T_warm;
static Real T_cut_mul;
static Real T_cut;

static Real g = 5./3.;


#endif
