import numpy as np

R = 1360.
n_fs = 1600
ndr = n_fs**0.5
dr = R/ndr

circum_t = 0.

for i in range(int(ndr)):
    circum_t += 2*np.pi * (dr+i*dr)

arc_l = circum_t/n_fs

c=0

for i in range(int(ndr)):
    circum = 2*np.pi * (dr+i*dr)
    ndtheta = circum/arc_l
    dtheta = 2*np.pi/ndtheta

    for j in range(int(ndtheta)):
        c+=1

print c
