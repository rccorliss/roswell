import numpy as np
import scipy.special as sp
from matplotlib import pyplot as plt

a = 20
b = 78
L = 105 

min_power = 5
max_power = 15

max_index = 0
number_zeros = 30

m = 0

start = 1.0e-6

zeros_list = []

while m <= max_index:
    n = 1
    zeros = []

    beta_list = []
    bessel_list = []

    beta = start

    def R_zero(beta):
        return sp.yv(m,a*beta)*sp.jv(m,b*beta) - sp.jv(m,a*beta)*sp.yv(m,b*beta)

    r_prev = R_zero(beta)

    if r_prev > 0:
        prev_sign = 1
    if r_prev < 0:
        prev_sign = -1

    while n <= number_zeros:
        beta += (10**(-min_power))
        r_new = R_zero(beta)
    
        if r_new > 0:
            new_sign = 1
        if r_new < 0:
            new_sign = -1
            
        if new_sign == prev_sign:
            zeros = zeros
        if new_sign != prev_sign:
            power = min_power + 1
            beta_sweep = beta - (10**(-min_power))
            r_prev_sweep = R_zero(beta_sweep)
                        
            while power <= max_power:
                new_sign_sweep = prev_sign
                step_size = 10**(-power)

                while new_sign_sweep == prev_sign:
                    beta_sweep += step_size
                    r_new_sweep = R_zero(beta_sweep)
    
                    if r_new_sweep > 0:
                        new_sign_sweep = 1
                    if r_new_sweep < 0:
                        new_sign_sweep = -1

                    if new_sign_sweep == prev_sign:
                        r_prev_sweep = r_new_sweep

                    if new_sign_sweep != prev_sign:
                        beta_sweep -= step_size
                        power += 1
            
            slope = (r_new_sweep - r_prev_sweep)/(10**(-max_power))
            root = beta_sweep - (r_prev_sweep/slope)
            zeros.append(root)
            
            n += 1
        
        prev_sign = new_sign
        r_prev = r_new

    print(m)
    print(zeros)
    zeros_list.append(zeros)

    m += 1

print(zeros_list)

#def Beta(m, n):
    #return zeros_array[m, n - 1]

#def Z(z, z_prime, m, n):
    #if z <= z_prime:
        #return (np.cosh(Beta(m, n) * z))*(np.sinh(Beta(m, n)*(L - z_prime)))/(np.sinh(Beta(m,n)*L))
    #if z > z_prime:
        #return -(np.cosh(Beta(m, n)* (L - z)))*(np.sinh(Beta(m, n)*(z_prime)))//(np.sinh(Beta(m,n)*L))
    
#def R(r, m, n):
    #return sp.yv(m,a*Beta(m, n))*sp.jv(m,r*Beta(m, n)) - sp.jv(m,a*Beta(m,n))*sp.yv(m,r*Beta(m, n))

#def N(m, n):
    #return (2/((((np.pi)**2))*((Beta(m, n))**2)))*((((sp.jv(m, a*Beta(m, n)))**2)/(((sp.jv(m, b*Beta(m, n)))**2))) - 1)

#def R_total(r, r_prime, m, n):
    #return (R(r, m, n)*R(r_prime, m, n))/N(m, n)

#def Phi(phi, phi_prime, m):
    #return np.cos(m * (phi - phi_prime))

#def Delta(m):
    #if m == 0:
        #return 1
    #else:
        #return 2

#def E_z(r, phi, z, r_prime, phi_prime, z_prime):

    #m = 0

    #E_z_array = np.zeros([max_index + 1, number_zeros])

    #while m <= max_index:
        #n = 1
        #while n <= number_zeros:
            #E_z_array[m, n - 1] = (1/(2*(np.pi)))*Delta(m)*Phi(phi, phi_prime, m)*R_total(r, r_prime, m, n)*Z(z, z_prime, m, n)
            #n += 1
        #m += 1

    #E_Z = np.sum(E_z_array)

    #return E_Z


