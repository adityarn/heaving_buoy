import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import linecache

c = 1
x = np.zeros(500)
y = np.zeros(500)
z = np.zeros(500)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for i in range(500):
    line = linecache.getline("xyz.txt",c)
    line_arr = line.split()

    x[c-1] = float(line_arr[0])
    y[c-1] = float(line_arr[1])
    z[c-1] = float(line_arr[2])
    
    c+=1

ax.scatter(x, y, z)
plt.show()
