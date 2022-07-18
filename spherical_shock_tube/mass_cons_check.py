from re import I
import numpy as np
import matplotlib.pyplot as plt


X_i = 100
X_f = 600

P_amb = 1
rho_amb = 1

gamma = 5/3

# Ambient sound speed
cs_amb = np.sqrt(gamma*P_amb/rho_amb)

# Incident shock front velocity
u_s = cs_amb/np.sqrt(gamma*X_i)

# Shocked medium sound speed
cs_shock = np.sqrt(gamma*P_amb/(X_f*rho_amb))

# Velocity of shocked medium
v = (1-X_i/X_f) * u_s

# Cloud of the radius
R = 0.5

t_arr = np.linspace(0,1.75,num=1000)

# Mass conservation LHS
Mcons_LHS = X_i*R**3

# Mass conservation RHS
Mcons_RHS  = X_i * (R - u_s*t_arr)**3
Mcons_RHS += X_f * ( (R - v*t_arr)**3 - (R - u_s*t_arr)**3  )

R_Mcons  = X_i*R**3/X_f
R_Mcons += (1-X_i/X_f) * (R-u_s*t_arr)**3
R_Mcons = R-R_Mcons**(1/3)

dR_Mcons = np.roll(R_Mcons,-1) - R_Mcons
dt = np.roll(t_arr,-1) - t_arr

v_Mcons = dR_Mcons/dt
v_Mcons = v_Mcons[:-1]


if __name__=="__main__":
    plt.figure()

    plt.plot(t_arr, Mcons_RHS)
    plt.axhline(Mcons_LHS, linestyle='dashed')

    plt.yscale('log') 

    plt.figure()

    plt.plot(t_arr,R-v*t_arr)
    plt.plot(t_arr,R-u_s*t_arr)

    plt.axhline(R,linestyle='dashed')
    plt.axhline(0,linestyle='dashed')


    plt.figure()
    plt.plot(t_arr[:-1], v_Mcons)
    plt.axhline(v,linestyle='dashed',color='tab:red')
    plt.axhline(0,linestyle='dashed',color='tab:red')

    # plt.yscale('log')
    # plt.ylim(1e-3,1e-1)
    plt.xlim(t_arr[0],t_arr[-1])

    plt.figure()
    plt.plot(t_arr, R_Mcons)
    plt.xlim(t_arr[0],t_arr[-1])