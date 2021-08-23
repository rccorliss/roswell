import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt

m = 20

def f1(x):
    return sp.jv(m, x + 1)

def f2(x):
    return sp.yv(m, x + 1)

x_min = 90.26354816250
x_max = 90.26354816251

y_min = f1(x_min)
y_max = f1(x_max)

slope = (y_max - y_min)/(x_max - x_min)

def y1(x):
    return slope * (x - x_min) + y_min

x = x_min
x_list = []
x1_list = []
x2_list = []

while x <= x_max:
    x1 = f1(x)
    #x2 = f2(x)

    x_list.append(x)
    x1_list.append(x1)
    #x2_list.append(x2)

    x += 1.0e-13


plt.plot(x_list, x1_list)
#plt.plot(x_list, x2_list)
plt.show()
    
