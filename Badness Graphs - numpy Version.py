import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt

np.seterr(divide='ignore')

min_power = 5
max_power = 15

number_zeros = 100

a = 20
b = 78
L = 105

zeros = []

start = 1.0e-6

def R_zero(beta, m):
    return sp.yv(m,a*beta) * sp.jv(m,b*beta) - sp.jv(m,a*beta) * sp.yv(m,b*beta)

m_list = []
delta_list = []

r = 65
r_prime = 30

m = 0

while m <= 50:
    
    zeros = []
    n = 0
    beta = start

    r_prev = R_zero(beta, m)

    if r_prev > 0:
        prev_sign = 1
    if r_prev < 0:
        prev_sign = -1

    while n <= number_zeros:
        beta += (10**(-min_power))
        r_new = R_zero(beta, m)
    
        if r_new > 0:
            new_sign = 1
        if r_new < 0:
            new_sign = -1
            
        if new_sign == prev_sign:
            zeros = zeros
        if new_sign != prev_sign:
            power = min_power + 1
            beta_sweep = beta - (10**(-min_power))
                        
            while power <= max_power:
                r_prev_sweep = R_zero(beta_sweep, m)

                if r_prev_sweep > 0:
                    prev_sign_sweep = 1
                if r_prev_sweep < 0:
                    prev_sign_sweep = -1

                new_sign_sweep = prev_sign_sweep

                while new_sign_sweep == prev_sign_sweep:
                    beta_sweep += (10**(-power))
                    r_new_sweep = R_zero(beta_sweep, m)
    
                    if r_new_sweep > 0:
                        new_sign_sweep = 1
                    if r_new_sweep < 0:
                        new_sign_sweep = -1
                    
                    if new_sign_sweep == prev_sign_sweep:
                        r_prev_sweep = r_new_sweep
                    if new_sign_sweep != prev_sign_sweep:
                        beta_sweep = beta_sweep - (10**(-power))
                        power += 1
            
            slope = (r_new_sweep - r_prev_sweep)/(10**(-max_power))
            root = beta_sweep - (r_prev_sweep/slope)
            zeros.append(root)

            print(n)
                                
            n += 1
        
        prev_sign = new_sign
        r_prev = r_new

    def Beta(n):
        return zeros[n - 1]

    def R(r, n):
        return sp.yv(m,a*Beta(n)) * sp.jv(m,r*Beta(n)) - sp.jv(m,a*Beta(n)) * sp.yv(m,r*Beta(n))
 
    def N(n):
        return (2/(((np.pi)**2)*((Beta(n))**2))) * ((((sp.jv(m, a*Beta(n)))**2)/(((sp.jv(m, b*Beta(n)))**2))) - 1)

    def R_total(r, r_prime, n):
        return (R(r, n)*R(r_prime, n))/N(n)

    def Delta_function(r, r_prime):
        count = 1
        nascent_delta = 0
        while count <= number_zeros:
            nascent_delta += R_total(r, r_prime, count)
            count += 1
        return nascent_delta
    
    delta = np.log(np.absolute(Delta_function(r, r_prime)))

    m_list.append(m)
    delta_list.append(delta)

    print(m)
    
    m += 1


plt.plot(m_list, delta_list)
plt.show()
