import numpy as np
import matplotlib.pyplot as plt

# M = 60
# N = 220
# K = 85 # target layer
M = 12
N = 44
K = 17

dx = 20
dy = 10
dz = 2

x = np.arange(M) * dx # for x coord
y = np.arange(N) * dy # for y coord

u = np.loadtxt('result.txt')
# u_reshaped = np.reshape(abs(u), (N, M))
u_reshaped = np.reshape(u, (N, M))

X, Y = np.meshgrid(x, y)
fig = plt.figure(figsize = (16, 8))
ax1 = fig.add_subplot(1, 2, 1)
cf = ax1.contourf(X, Y, u_reshaped, 100, cmap = 'gnuplot2')
fig.colorbar(cf, ax = ax1)

fig.show()