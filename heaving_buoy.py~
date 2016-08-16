import numpy as np
import linecache


'''Reading From Input File "xyz.txt"
which contains all coordinates, normals
and area of each panel. This file is
generated with the C++ code "coord_gen.cpp" '''

c = 0
line = linecache.getline("xyz.txt",c+1)
line_arr = line.split()

Time_period = float(line_arr[0])
L = float(line_arr[1])
n_total = int(line_arr[2])
n_body = int(line_arr[3])
n_fs = int(line_arr[4])
n_bottom = int(line_arr[5])
n_ff = int(line_arr[6])

Kl = 2*np.pi/L
sigma = 2*np.pi/Time_period

x = np.zeros(n_total)
y = np.zeros(n_total)
z = np.zeros(n_total)
nx = np.zeros(n_total)
ny = np.zeros(n_total)
nz = np.zeros(n_total)
ar = np.zeros(n_total)
a = np.zeros((n_total,n_total),dtype=complex)
b = np.zeros(n_total,dtype=complex)
phi = np.zeros(n_total,dtype=complex)

for c in range(n_total):
    line = linecache.getline("xyz.txt",c+2)
    line_arr = line.split()
    x[c] = float(line_arr[0])
    y[c] = float(line_arr[1])
    z[c] = float(line_arr[2])
    nx[c] = float(line_arr[3])
    ny[c] = float(line_arr[4])
    nz[c] = float(line_arr[5])
    ar[c] = float(line_arr[6])

'''Green's mixed Distribution coefficients and RHS matrix are 
computed here'''
for i in range(n_total):
    for k in range(n_body):
        if(i!=k):
            r = (x[i]-x[k])**2 + (y[i]-y[k])**2 + (z[i]-z[k])**2
            numer = -sigma * nz[k]* ar[k]
            denom = r**0.5
            b[i] += complex(0,numer/denom)

    for k in range(n_body):
        if(i != k):
            r = (x[i]-x[k])**2 + (y[i]-y[k])**2 + (z[i]-z[k])**2
            numer = ((x[i]-x[k])*nx[k] + (y[i]-y[k])*ny[k] + (z[i]-z[k])*nz[k])*ar[k]
            denom = r**1.5
            a[i][k] = complex(numer/denom,0)
        else:
            a[i][k] = complex(2*np.pi,0)
    
    for k in range(n_body,n_body+n_fs,1):
        if(i != k):
            r = (x[i]-x[k])**2 + (y[i]-y[k])**2 + (z[i]-z[k])**2
            numer = (x[i]-x[k])*nx[k] + (y[i]-y[k])*ny[k] + (z[i]-z[k])*nz[k]
            denom = r**1.5
            numer2 = sigma**2
            denom2 = 9.81 * r**0.5
            a[i][k] = complex((numer/denom - numer2/denom2)*ar[k],0)
        else:
            a[i][k] = complex(2*np.pi,0)

    for k in range(n_body+n_fs,n_body+n_fs+n_bottom,1):
        if(i != k):
            r = (x[i]-x[k])**2 + (y[i]-y[k])**2 + (z[i]-z[k])**2
            numer = ((x[i]-x[k])*nx[k] + (y[i]-y[k])*ny[k] + (z[i]-z[k])*nz[k])*ar[k]
            denom = r**1.5
            a[i][k] = complex(numer/denom,0)
        else:
            a[i][k] = complex(2*np.pi,0)
            
    for k in range(n_body+n_fs+n_bottom,n_total,1):
        if(i != k):
            r = (x[i]-x[k])**2 + (y[i]-y[k])**2 + (z[i]-z[k])**2
            numer = ((x[i]-x[k])*nx[k] + (y[i]-y[k])*ny[k] + (z[i]-z[k])*nz[k])*ar[k]
            denom = r**1.5
            numer2 = Kl/r
            a[i][k] = complex((numer/denom*ar[k]), -numer2*ar[k])
        else:
            a[i][k] = 2*np.pi

phi = np.linalg.solve(a,b)

fjk = complex(0.0,0.0)
for i in range(n_body):
    fjk +=  phi[i] * nz[i] * ar[i]

fjk = fjk * 1025. * sigma * complex(0,-1)

mu = fjk.real /(-sigma**2)
damping = fjk.imag/(-sigma)

'''Output the results into a text file'''
print"\n\n***********************************************************\n\nGreen's Therem solved for Heaving Buoy\n"
print "n_total = %d" %(n_total) + "\nn_body = %d"%n_body + "\nn_bottom = %d"%n_bottom + "\nn_fs = %d"%n_fs +"\nn_ff = %d"%n_ff
print "\nfjk = %s"%fjk + "\nmu = %f"%mu + "\ndamping = %f"%damping

f = open("run_info",'a')
f.write("\n\n***********************************************************\n\nGreen's Therem solved for Heaving Buoy\n")
f.write("n_total = %d" %(n_total) + "\nn_body = %d"%n_body + "\nn_bottom = %d"%n_bottom + "\nn_fs = %d"%n_fs +"\nn_ff = %d"%n_ff)
f.write("\nfjk = %s"%fjk + "\nmu = %f"%mu + "\ndamping = %f"%damping)
f.close()
