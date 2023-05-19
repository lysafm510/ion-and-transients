import numpy as np
import linecache
import os
import csv

# path = "20220501_random4"

data = np.loadtxt("avg_gn_jsr_gauss_new_2.csv")
# return len(data) - 3
# return len(data) - 2

a = data[0]
n = data.size
with open("avg_gn_jsr_new_eee.csv", "w", newline='') as file_object:
    writer = csv.writer(file_object)
    for i in range(n):
        b = data[i] / a
        writer.writerow([b])
