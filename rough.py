import numpy as np

n_body = 1800

ndy = int((n_body*.5)**0.5)
ndtheta = ndy
dtheta = 2*np.pi/ndtheta
dy = 3.0/ndy
