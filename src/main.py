
import numpy as np
import scipy.linalg as sla
from scipy import sparse 
import matplotlib.pyplot as plt
import seaborn as sns
import time


def parse1(filename, M, N, K, Kc):
    Dx  = []
    Dy = []
    # i = 0
    first_layer_x = (M * N * Kc)/6
    last_layer_x = (M * N * (Kc + 1))/6
    first_layer_y = (M * N * (K + Kc))/6 
    last_layer_y = (M * N * (K + Kc + 1))/6
    
    with open(filename) as file:
        # print("file was onened")
        data = file.read().splitlines()
        for id, line in  enumerate (data):
            row = line.split('\t')
            for elem in row:
                value = elem.strip()
                if(id >= first_layer_x  and  id < last_layer_x):
                    if value == '':
                        continue
                    else:
                        Dx.append(float(value))
                if(id > first_layer_y  and  id <= last_layer_y):
                    if value == '':
                        continue
                    else:
                        Dy.append(float(value))                       
    return(Dx, Dy)


def k_right(i):
    return 2 * Dx[i] * Dx[i + 1] / (dx ** 2 * (Dx[i] + Dx[i + 1]))


def k_left(i):
    return 2 * Dx[i] * Dx[i - 1] / (dx ** 2 * (Dx[i] + Dx[i - 1]))


def k_up(i):
    return 2 * Dy[i] * Dy[i - M] / (dy ** 2 * (Dy[i] + Dy[i - M]))


def k_down(i):
    return 2 * Dy[i] * Dy[i + M] / (dy ** 2 * (Dy[i] + Dy[i + M]))


def solve(M, N, K, Dx, Dy):

    R = [0] * (M * N)
    A = [[0] * (M * N) for i in range(M * N)]

    counter_left = 0
    counter_right = M - 1

    for i in range(M * N):
        if i < M:
            R[i] += 2 * u_y0 * (Dy[i] / (dy ** 2))
            A[i][i] += 2 * Dy[i] / (dy ** 2)
        else:
            R[i] += 0
            A[i][i - M] -= k_up(i)
            A[i][i] += k_up(i)

        if i > M * (N - 1) - 1:
            R[i] += 2 * u_yN * (Dy[i] / (dy ** 2))
            A[i][i] += (2 * (Dy[i] / (dy ** 2)))

        else:
            R[i] += 0
            A[i][i + M] -= k_down(i)
            A[i][i] += k_down(i)

        if i == counter_left:
            R[i] += 2 * u_x0 * (Dx[i] / (dx ** 2))
            A[i][i] += (2 * (Dx[i] / (dx ** 2)))
            counter_left += M
        else:
            R[i] += 0
            A[i][i - 1] -= k_left(i)
            A[i][i] += k_left(i)

        if i == counter_right:
            R[i] += 2 * u_xM * (Dx[i] / (dx ** 2))
            A[i][i] += (2 * (Dx[i] / (dx ** 2)))
            counter_right += M
        else:
            R[i] += 0
            A[i][i + 1] -= k_right(i)
            A[i][i] += k_right(i)
    return A, R


if __name__ == '__main__':
    
#   params of mesh:
    M = 60
    N = 220
    K = 85
#   target layer
    Kc = 10
#   boundary dirichlet
    u_x0 = 0
    u_xM = 10
    u_y0 = 1
    u_yN = 15
#   grid step (size of cells)
    dx = 20
    dy = 10
    dz = 2
    
    # start = time.time()
    Dx,Dy = parse1("data/spe_perm.dat", M, N, K, Kc)
    # print('parsing was done')
    # end = time.time() - start
    # print("Time LENA", end)

    A, R = solve(M, N, K, Dx, Dy)
    # print('A-matrix and b-vector are made')

#   matrix drawing
    # plt.figure(figsize=(10, 10))
    # sns.heatmap(A)
    # plt.show()
    
#   various solver options
    res = np.linalg.solve(A, R)
    # res = sla.solve(A, R)
    
    f = open("result.txt", "w")
    for val in res:
        f.write(str(val) + "\n")
    f.close()   