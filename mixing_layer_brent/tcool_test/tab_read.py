import numpy as np
import matplotlib.pyplot as plt

plt.style.use('../../plot_scripts/plot_style.mplstyle')

# import globals as g

CONST_pc  = 3.086e18
CONST_yr  = 3.154e7
CONST_amu = 1.66053886e-24
CONST_kB  = 1.3806505e-16
unit_length = CONST_pc*1e3  # 1 kpc
unit_time   = CONST_yr*1e6  # 1 Myr
unit_density = CONST_amu    # 1 mp/cm-3
unit_velocity = unit_length/unit_time
KELVIN = unit_velocity*unit_velocity*CONST_amu/CONST_kB

Xsol = 1.0
Zsol = 1.0

X = Xsol * 0.7381
Z = Zsol * 0.0134
Y = 1 - X - Z

mu  = 1.0/(2.*X+ 3.*(1.-X-Z)/4.+ Z/2.);
mue = 2.0/(1.0+X);
muH = 1.0/X;
mH = 1.0


x_ind   = 1
rho_ind = 2
prs_ind = 3
vel_ind = 4

gamma = 5/3
T_floor = 4e4
T_hot   = 4e6

cool_arr  = np.loadtxt('cool_coef.txt')
cool_t    = cool_arr[:,0]
cool_coef = cool_arr[:,1]

data1 = np.loadtxt(f'Turb.block0.out2.00000.tab')
data2 = np.loadtxt(f'Turb.block0.out2.00001.tab')

# dt = 1.135926e-03 * 5
dt = 1e-3 
Lambda0 = 2

x1   = data1[:,x_ind]
rho1 = data1[:,rho_ind]
prs1 = data1[:,prs_ind]
vel1 = data1[:,vel_ind]

x2   = data2[:,x_ind]
rho2 = data2[:,rho_ind]
prs2 = data2[:,prs_ind]
vel2 = data2[:,vel_ind]

T1 = (prs1/rho1)*KELVIN*mu
T2 = (prs2/rho2)*KELVIN*mu

TE1 = prs1/(gamma-1) * unit_density * (unit_velocity**2)
TE2 = prs2/(gamma-1) * unit_density * (unit_velocity**2)

dTE = (TE2-TE1)/ (dt * unit_time) 
Lambda = -dTE/(rho1**2) 
dT = T2-T1

#*_________________________________

plt.figure(figsize=(10,10))

plt.plot(T1,Lambda/1e-23, label='Simulation' )
# plt.plot(T1,x1, label='Simulation' )
plt.plot(cool_t,cool_coef*Lambda0, label='From table')

plt.xlim(T_floor, T_hot)

plt.title(r'$\Lambda_0$ = 1')

plt.xlabel('T (K)')
plt.ylabel(r'$\Lambda(\rm T)$ (10$^{-23}$ ergs cm$^{3}$/s)')

plt.axvline(0.5*T_hot, linestyle='dotted', label=r'T$_{\rm cut}$')

plt.legend()

plt.xscale('log')
plt.yscale('log')
#*_________________________________

plt.figure(figsize=(10,10))

plt.scatter(x1,T2, label='P2' )
plt.scatter(x1,T1, label='P1' )
# plt.plot(x2,T2, label='Final temperature' )

plt.axhline(T_floor, linestyle='dotted')
plt.axhline(T_hot, linestyle='dotted')

plt.title(r'$\Lambda_0$ = 1')

plt.xlabel('x')
plt.ylabel('P')

plt.legend()

# plt.xscale('log')
# plt.yscale('symlog')
plt.yscale('log')
#*_________________________________

plt.figure(figsize=(10,10))

plt.scatter(x1,vel2, label='P2' )

plt.title(r'$\Lambda_0$ = 1')

plt.xlabel('x')
plt.ylabel('P')

plt.legend()

# plt.xscale('log')
# plt.yscale('symlog')
# plt.yscale('log')
#*_________________________________