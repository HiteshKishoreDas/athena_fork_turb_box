from math import gamma
import numpy as np

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

# Reflected shock front velocity
u_r = -v + np.sqrt(v**2 + 4*cs_amb**2/(gamma*X_f))
u_r /= 2

# Shock2-ed medium sound speed
rho_1 = P_amb/u_r**2
cs_shock_2 = np.sqrt(gamma*P_amb/(rho_1))


print(f'cs_shock  : {cs_shock}')
print(f'u_s       : {u_s}')
print(f'u_r       : {u_r}')
print(f'cs_shock_2: {cs_shock_2}')