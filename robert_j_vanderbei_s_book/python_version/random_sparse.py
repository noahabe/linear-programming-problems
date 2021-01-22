from numpy import *
import numpy as np
from utils import *
import pd 
import twophase
import hsd
import hsdls
import intpt

m = 40
n = 30

a = zeros(m*n)
b = zeros(n+m)
c = zeros(n+m)
ia = zeros(m*n).astype(int16)
ka = zeros(n+m+1).astype(int16)
w = zeros(m)
x = zeros(n)
y = zeros(m)
z = zeros(n)
Ax = zeros(m)
Aty = zeros(n)

k = 0
for j in range(0,n):
    ka[j] = k
    colden = np.random.uniform()
    colden = colden**2
    ii = int(np.random.uniform()*m)
    ia[k] = ii
    a[k]  = round(20*(np.random.uniform()-0.5))
    k = k+1
    for i in range(0,m):
        if (i != ii):
            if (np.random.uniform() < colden):
                ia[k] = i
                a[k]  = round(20*(np.random.uniform()-0.5))
                k = k+1
nz = k
ka[n] = k

for j in range(0,n):
    x[j] = round(10*np.random.uniform())
    z[j] = round(10*np.random.uniform())

for i in range(0,m):
    w[i] = round(10*np.random.uniform())
    y[i] = round(10*np.random.uniform())

for i in range(0,m): Ax[i]  = 0
for j in range(0,n): Aty[j] = 0
for j in range(0,n):
    for k in range(ka[j],ka[j+1]):
        i = ia[k]
        Ax[i]  = Ax[i]  + a[k]*x[j]
        Aty[j] = Aty[j] + a[k]*y[i]

for i in range(0,m): b[i] = Ax[i]  + w[i]
for j in range(0,n): c[j] = Aty[j] - z[j]

f = 0

x = zeros(n+m)
y = zeros(n+m)
w = zeros(n+m)
z = zeros(n+m)

for j in range(0,10):
    print('j = ',j)
    print("Primal Dual Simplex (pd)")
    status = pd.solver(m, n, nz, ia, ka, a, b, c, f, x, y, w, z)
    print(statmsg[status])
    pd.pd_reset()

print("Two Phase (2phase)")
status = twophase.solver(m, n, nz, ia, ka, a, b, c, f, x, y, w, z)
print(statmsg[status])
print()
pd.pd_reset()

print("Homegeneous Self-Dual (hsd)")
status = hsd.solver(m, n, nz, ia, ka, a, b, c, f, x, y, w, z)
print(statmsg[status])
print()
pd.pd_reset()

print("Homegeneous Self-Dual LS (hsdls)")
status = hsdls.solver(m, n, nz, ia, ka, a, b, c, f, x, y, w, z)
print(statmsg[status])
print()
pd.pd_reset()

print("Interior-Point (intpt)")
status = intpt.solver(m, n, nz, ia, ka, a, b, c, f, x, y, w, z)
print(statmsg[status])
print()
pd.pd_reset()
