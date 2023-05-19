import random
from math import sqrt, log

import numpy as np

n = []
for i in range(100):
    h = random.gauss(0, 1)
    n.append(h)

arr = np.array(n)

u = sqrt(-2 * log(1 / 3))
arr = arr[(arr > -u) & (arr < u)]
arr = arr * 100 / u + 100
arr = arr.astype(int)
arr = np.sort(arr)
print(arr)
print(arr.std())

# print(arr.std())

m = []
for i in range(100):
    o = random.gauss(0, 1)
    m.append(int(o * 100 / u) + 100)

m = np.array(m)
m = np.sort(m)
m = m[(m > 0) & (m < 200)]
print(m)
print(m.std())
