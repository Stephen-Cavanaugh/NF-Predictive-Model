# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 16:08:11 2022
@author: Stephen
"""
import pandas as pd
import os
import numpy as np
path = os.getcwd()

import numpy as np
import scipy as sp
import matplotlib as mp
import matplotlib.pyplot as plt
from scipy.integrate import quad
import sympy as syp
from sympy.solvers import solve
from sympy import Symbol, exp
from scipy.optimize import minimize_scalar
import math
import sys


lit_reject = pd.read_excel('C:/your_path_/your_file.xlsx', sheet_name = 'NTR_7250')

compound_name = np.array([x for x in lit_reject['Compound']])
membrane = np.array([x for x in lit_reject['membrane']])

lit_reject["comp_name_mem"] = lit_reject["Compound"] + "_" + lit_reject["membrane"]

comp_name_mem = np.array([x for x in lit_reject['comp_name_mem']])

charge = np.array([x for x in lit_reject['charge']])

pore_sz_lit = np.array([x for x in lit_reject['pore size (nm)']])

pore_sz_plot = pore_sz_lit*10**(-9)

mem_charge_lit = np.array([x for x in lit_reject['mem charge']])

length_lit = np.array([x*(10**(-9)) for x in lit_reject['Molar length']])

width_lit = np.array([x*(10**(-9)) for x in lit_reject['molar width']])

molar_mass = np.array([x for x in lit_reject['molar weight']])

conc_lit = np.array([x for x in lit_reject['concentration']])

pole_lit = np.array([x for x in lit_reject['zwit']])

Rej_obs = np.array([x for x in lit_reject['rejection']])

N_steps = 9

Final_rej_array = np.zeros(len(length_lit))

Rej_pore_array = [[None for y in range(len(lit_reject))] for x in range(len(lit_reject)+1)]

for m in range(0,len(compound_name)):
    sigma_star_array = 0.15E-9 #standard deviation of the normal dist
    r_min_array = (pore_sz_lit[m]*10**(-9))
    r_max_array = 2.5*(pore_sz_lit[m]*10**(-9)) #maximum pore size (to truncate the distribution) - m
    Pore_array_steps = np.linspace(r_min_array,r_max_array,N_steps)
    
    
    length = length_lit[m] #Length of the compound - from Kiso et al 2010. Will use chemsketch or other to determine width of each compound as necessary - m
    width = width_lit[m] #width of the compound - from Kiso et al 2010. Will use chemsketch or other to determine width of each compound as necessary - m
    
    pressure = 72.5*6895 #Pa 
    approach = np.linspace(0,2*np.pi,N_steps) #approach angle, in radians
    R_gas = 8.314 #Ideal gas constant - j/molK
    T = 293.13 #Temp - K
    F = 96485.3 #Faraday constant - C/mol
    M_weight = molar_mass[m] #Molecular weight of the compound of interest - g/mol
    C_mass_mg = conc_lit[m]#Concentration in the bulk in mg/L
    Cmass = C_mass_mg/1000/M_weight #Concentration in the bulk in mol/L

    K = 1.380649 * 10**(-23) #Boltzman constant - J/K

    """
    stokes radius, which will be used later for the 'intra-pore' hindrance calcs
    """
    r_stokes = (1.42*width/(1E-9)-0.142)*1E-9

    r_stokes_diff = (1.42*width/(1E-9)-0.142)*1E-9
    d_water_diff = 0.28E-9
    eta_water = 0.001002
    
    thickness = 6.4E-8 #Thickness of the membrane - m #From "Development of Anti-fouling membranes for water treatment" - Steven Weinman 2018

    C_bulk = (Cmass)*1000 #Concentration in the bulk - mol/m3

    e = 1.60217662*10**(-19) #elementary charge - Coulomb

    """
    Defining some arrays that will be used later for plotting. array of phi, rejection, and aproach (alpha)
    """ 
    phi_array = np.zeros(N_steps)
    
    
    Rej_pore_inner = [None for x in range(N_steps)]
    alpha_array = np.linspace(0,np.pi/2,N_steps)

    pH = 3

    
    zeta = mem_charge_lit[m]


    
    Z_array = np.array([])

    """
    for loop that will step through all possible orientation that the compound can approach from. 
    may need to edit slightly in a future code so that this is a distribution? as is, it assumes
    all are oriented the same way. may be able to average rejections from all orientations if we 
    assume it is equally likely that it approaches from any given angle? 
    """
    pole = pole_lit[m]
    
    pole_charge = charge[m]
    
    pole_charge_zwit = charge[m]
    
    for p in range(0,len(Pore_array_steps)):
        R_array = []
        eta_diff = eta_water*(1+18*(d_water_diff/Pore_array_steps[p])-9*(d_water_diff/Pore_array_steps[p])**2)
        Diff = (K*T/(6*np.pi*eta_diff))*(1/r_stokes_diff)#1E-9 #diffusivity - placeholder m2/s
        
        pore_sz = Pore_array_steps[p] #0.69E-9 #m # see 10.1002/9780470024539.ch4 fig 5 for conversion from kda to nm
        
        for i in range(len(approach)):
            if pole == 0:
                N_ion = 0
                z_i = 0
            elif pole == 1:
                if 0 <= approach[i] < np.pi/2 or np.pi <= approach[i] < np.pi*2:
                    z_i = pole_charge
                    N_ion = pole_charge
                else:
                    z_i = 0
                    N_ion = 0
            else:
                if 0 <= approach[i] <= np.pi:
                    z_i = pole_charge
                    N_ion = pole_charge
                else:
                    z_i = pole_charge_zwit
                    N_ion = pole_charge_zwit
            Z_array = np.append(Z_array, z_i)
                        
            """
            Defines a function to integrate to find the final phi term. From
            Kiso et al  
            """
            def integrand_phi(alpha,phi_alpha):
                return np.sin(alpha)*phi_alpha
            
            L = length 
            MW = width/2 
            r_pore_array = pore_sz
    
            L_pore = L*np.cos(approach[i])+MW*np.sin(approach[i])
            a_width = np.sqrt(r_pore_array**2 - ((L_pore/2)**2)) - MW
            phi_alpha = (a_width/r_pore_array)**2
            
            """
            'quad' calls the integration of a function - in this case 'integrand_phi' - from 0 to np.pi/2, with 
            function arguments l, r_pore, and mw.
            """
            phi = quad(integrand_phi, 0,(np.pi/2), args=(phi_alpha))#L,r_pore_array,MW))
            
                    
            phi_array[i] = phi[0] #defining our final phi term based on the integration above. 'quad' returns an array, so we need to pull the first value (0 in python) from that array
            
            r_cav = r_stokes*1.07 #Rahin & Honig 1985, Reevaluation of hte Born model of ion hydration, and Szymczyk 2007 
            
            """
            setting up the hindrance factors
            """ 
            lamb = r_stokes/r_pore_array #Equation 5
            
            G_lamb = 1+0.054*lamb -0.988*lamb**2+0.441*lamb**3 #euqation 9
            K_inv = 1-2.3*lamb+1.154*lamb**2+0.224*lamb**3 #equation 8
            k_id = K_inv #equation 6
            
            
            phi_old = phi_array[i]
            k_ic = (2-phi_old)*G_lamb #equation 7
            
            
            """
            Defining some constants- primarily found from Bowen and Wellfoot 2002, 
            Modelling of membrane nanofiltration - pore size distribution effects
            """
            eta_water = 0.001002 #bulk viscosity of water - N s/m^2
            
            r_star = pore_sz #average pore size - m
            sigma_star = 0.15E-9 #standard deviation of the normal dist
            d_water = 0.28E-9 #diameter of a water molecule - nm
            r_max = 1.25*r_star #maximum pore size (to truncate the distribution) - m
            
            
            """
            average viscosity increase due to the layer of water inside the pore, Bowen and wellfoot 2002,
            Modelling of membrane nanofiltration - pore size distribution effects, equation 3
            """
            eta_eta0 = 1+18*(d_water/r_pore_array)-9*(d_water/r_pore_array)**2 #average viscosity increase due to the layer of water inside the pore
            
            
            """
            Diffusion coefficient, adjusted for hindrance effects, Kiso et al, 2010, equation 2 - m2/s
            """
            D_inf = k_id*Diff*(eta_eta0**(-1)) 
            
            
            """
            The following two functions define the Jv and Ak terms, found in Bowen and Welfoot, 2002, Modelling of 
            membrane nanofiltration
            """
            def integrand_Jv(r,sigma_star,r_star,eta_water,d_water):
                return ((1/(r*np.sqrt(2*np.pi*np.log(1+(sigma_star/r_star)**2))))*(np.exp(-((np.log(r/r_star)+((np.log(1+(sigma_star/r_star)**2))/2))**2)/(2*(np.log(1+(sigma_star/r_star)**2)))))*(r**4))/(eta_water*(1 + (18*(d_water/r)) - 9*(d_water/r)**2))
            
            
            def integrand_ak(r,r_star,sigma_star):
                return (1/(r*np.sqrt(2*np.pi*np.log(1+(sigma_star/r_star)**2))))*(np.exp(-((np.log(r/r_star)+((np.log(1+(sigma_star/r_star)**2))/2))**2)/(2*(np.log(1+(sigma_star/r_star)**2)))))*r**2
            
            
          
            if z_i != 0:
                def f (x,c_i_rung): #defining our primary equation - eq 1
                    return (V_r/(k_id*D_inf))*(k_ic*c_i_rung - c_perm[q]) - ((z_i*F*c_i_rung)/(R_gas*T))*((((z_i*V_r)/(k_id*D_inf))*(k_ic*c_i_rung-c_perm[q]))/((F/(R_gas*T))*(c_i_rung*z_i**2)))
            else:
                def f (x,c_i_rung): #If there's no ionic charge, we remove the charge portion of the equation
                    return (V_r/(k_id*D_inf))*(k_ic*c_i_rung - c_perm[q])
            
                
            """
            Runge kutta fehlberg, following algorithm in Numerical Analysis, 9th, Burden and Faires, chapter 5. 
            This defines the function, which will be called shortly
            """
            def rung_kut_fehl4 (a, b, x_0, tol, hmax, hmin): 
                # global w
                # global h
                # global T
                global X
            
                FLAG = 1
                t = a #Distance along membrane pore
                w = x_0 #Concentration 
                h = hmax
                T1 = np.array([t])
                X = np.array([w])
                
                    
                while FLAG == 1:
                    """
                    Defining all of coefficients used later to make it a bit easier to construct. 
                    """    
                    K2_1 = 0.25
                
                    K3_1 = 3/8
                    K3_2 = 3/32
                    K3_3 = 9/32
                    
                    K4_1 = 12/13
                    K4_2 = 1932/2197
                    K4_3 = 7200/2197
                    K4_4 = 7296/2197
                    
                    K5_1 = 439/216
                    K5_2 = 3680/513
                    K5_3 = 845/4104
                    
                    K6_1 = 1/2
                    K6_2 = 8/27
                    K6_3 = 3544/2565
                    K6_4 = 1859/4104
                    K6_5 = 11/40
                    
                    R1 = 1/360
                    R2 = 128/4275
                    R3 = 2197/75240
                    R4 = 1/50
                    R5 = 2/55
                    
                    W1 = 25/216
                    W2 = 1408/2565
                    W3 = 2197/4104
                    W4 = 1/5
                    # global K1
                    # global K2
                    # global K3
                    # global K4
                    # global K5
                    # global K6
                    
                    """
                    Defining the individual values that get plugged into the R (stepping) term. each adds a small 'step' to the function
                    """
                    K1 = h*f(t,w)
                    K2 = h*f((t + K2_1*h), (w + K2_1*K1))
                    K3 = h*f((t + K3_1*h), (w + K3_2*K1 + K3_3*K2))
                    K4 = h*f((t + K4_1*h), (w + K4_2*K1 - K4_3*K2 + K4_4*K3))
                    K5 = h*f((t + h), (w + K5_1*K1 - 8*K2 + K5_2*K3 - K5_3*K4))
                    K6 = h*f((t + K6_1*h), (w - K6_2*K1 + 2*K2 - K6_3*K3 + K6_4*K4 - K6_5*K5))
                    
                    
                    R = 1/h *(np.abs(R1*K1 - R2*K3 - R3*K4 + R4*K5 + R5*K6))
                    print(R)
                    
                    """
                    Does some checking of the function
                    """
                    
                    if R <= tol:
                        t = t + h
                        w = w + W1*K1 + W2*K3 + W3*K4 - W4*K5
                        T1= np.append(T1, t)
                        X = np.append(X, [w], axis = 0)
                           
                        #return [t,w,h]
                    
                    delta = 0.84*(tol/R)**(1/4)
                    
                    
                    if delta <= 0.1:
                        h = 0.1*h
                        
                    elif delta >=4:
                        h = 4*h
                        
                    else:
                        h = delta*h
                        
                        
                    if h> hmax:
                        h = hmax
                    
                    if t>= b:
                        FLAG = 0
                        
                    elif (t+h) > b:
                        h = b - t
                        
                    elif h < hmin:
                        FLAG = 0
                        print('minimum h exceeded')
                        break     
                    
                return (T1, X)
                
            
            """
            Runge kutta of order 4, following algorithm in Numerical Analysis, 9th, Burden and Faires, chapter 5 (not used)
            """
            def rung_kut4 (x0, ci0, xn, n): 
                step = (xn - x0)/n
                print('\n--------SOLUTION--------')
                print('-------------------------')    
                print('x0\ty0\tyn')
                print('-------------------------')
                global x_steps 
                x_steps = np.zeros(n)
                global c_steps 
                c_steps = np.zeros(n)
                c_steps[0] = ci0
                
                for i in range(n):
                    
                    k1 = step * (f(x0, ci0))
                    k2 = step * (f((x0+step/2), (ci0 + k1/2)))
                    k3 = step * (f((x0+step/2), (ci0 + k2/2)))
                    k4 = step * (f((x0+step), (ci0 + k3)))
            
                    k = (k1 + 2*k2 + 2*k3 + k4)/6
                    
                    cn = ci0 + k 
                    print('%f\t%f\t%f'% (x0,ci0,cn) )
                    print('-------------------------')
                    ci0 = cn
                    x0 = x0+step
            
                    x_steps[i] = x0
                    c_steps[i] = cn
                    #print(x0,ci0,cn)
            
                return cn    
            
                print('\nAt x=%f, y=%f' %(xn,cn))
            
            x0 = 0 #starting point of the membrane - 0m through the pore essentially 
            xn = thickness #ending point of the membrane - m thickness. Used later on to define the bounds of the runge kutta
            
            """
            Somre more constants defined. 
            See bowen and welfoot 2002 - "modeling the performance of membrane nanofiltration - 
            critical assessment and model development" for epsi_b, and equation for epsi_p
            epsi_m is also defined here, noted as epsi_star in this paper
            epsi_b is assumed to be 80 (Szymczyk and Fievet 2005, and above bowen). 
            any range from 60-80 seems common. 
            """
            epsi_0 = 8.85419*10**(-12) #the vacuum permittivity - Physics of Chemistry of Interfaces - Butt and Kappl, 2003, Wiley
            
            epsi_m = 6 #Pulled from Szymczyk and Fievet 2005 - Investigating transport properties of nanofiltration membranes by means of a steric, electric and dielectric exclusion model. 
            
            epsi_b = 80 
            
            z_1 = z_i
            
            epsi_p = epsi_b - 2*(epsi_b-epsi_m)*(d_water/r_pore_array)+(epsi_b-epsi_m)*(d_water/r_pore_array)**2 #not in cheat sheet - see Bowen and Wellfoot 2002 -"modeling the performance of membrane nanofiltration - critical assessment and model development"
                
            del_wi = (((z_1*e)**2)/(8*np.pi*epsi_0*K*T*r_cav))*((1/epsi_p)-(1/epsi_b)) #equation 21
            
          
            C_bulk_avog = C_bulk*6.0221409E23
            
            del_psi = zeta - zeta * np.exp(-np.sqrt((2*C_bulk_avog*e**2)/(epsi_0*epsi_p*K*T))*d_water) #poisson boltzman equation, see "Pysics and chemistry of interfaces" Butt, Graf, and Kappl - equation 4.8 and 4.9
                
            
            """
            concentration just inside the pore, as defined by equation 20 in the cheat sheet. 
            """
            c_i_20 = C_bulk * phi_old*np.exp(-(F*z_i*del_psi + del_wi)/(R_gas*T))
            
            
          
            last_cperm = None
            
            q = 0
            c_perm  =  np.zeros(1)
    
            del_c_app = np.array([])
            
           
            while True:
                
                               
                """
                calling the integration of Jv from earlier using 'quad' function. Integrating from 0 to r_max, over r, and inputs the args as needed
                """
                
                
                """
                calling the integration of Ak from earlier using 'quad' function. Integrating from 0 to r_max, over r, and inputs the args as needed
                """
               
                
                """
                solving for the change in concentration (used in eq 13, which is then used in eq 11, which finally is used in Jv_full below) 
                (it is not included in the quad integration above, since it is effectively a 'constant' - at least as far as the integration is concerned
                 - it has no 'r' terms in it.) (it is also not used in the above function since it is dependent on concentration, which needs to be individually
                calculated during this loop)
                """
                del_c = C_bulk - c_perm[q]#c_i_20 - c_perm[q]
                
                del_c_app = np.append(del_c_app,del_c)
                
                R_g_pi = 8.31E3 # Pa/Kmol/L
                del_c_pi = C_mass_mg #mg/L
                M_weight_pi = M_weight*1000 #mg/mol
                
                del_pi = (N_ion*R_g_pi*del_c*T)/M_weight_pi
            
                del_pe = pressure - del_pi #equation 12
                
                V_r = (r_pore_array**2)*del_pe*(1/thickness)*(1/eta_diff)
                
                """
                finally calling the runge kutta fehlberg algorithm from above, using the arguments:
                x0 - bulk side of the pore
                xn - permeate side of the pore
                c_i_20 - initial concentration at the pore interface
                the tolerance allowed
                max and min step sizes allowed
                """
                rung = rung_kut_fehl4(x0,xn,c_i_20,(10**(-20)),xn/500,xn/10000)
            
                rung
                
                iter_rung = len(rung[1])
                
                """
                iterating on the concentration, and appending the last value to the 
                """
                c_perm = np.append(c_perm, rung[1][iter_rung-1])
            
                print('\nAfter loop %f, the value of permeate is: %f'% (q+1,  c_perm[q]))
            
                if last_cperm is not None and math.isclose(c_perm[q-1], last_cperm):
                    break
            
                q += 1 
            
                last_cperm = c_perm[q]
    
    
            Rej = (1- (c_perm[q-1]/C_bulk))*100 #standard rejection equation, not shown in cheat sheet
            R_array = np.append(R_array,Rej)
        
        Rej_pore_inner = np.vstack((Rej_pore_inner,R_array))

    print(R_array)

    """
    some boilerplate code to make the plot show pi in the x-axis
    """
    def multiple_formatter(denominator=2, number=np.pi, latex='\pi'):
        def gcd(a, b):
            while b:
                a, b = b, a%b
            return a
        def _multiple_formatter(x, pos):
            den = denominator
            num = int(np.rint(den*x/number))
            com = gcd(num,den)
            (num,den) = (int(num/com),int(den/com))
            if den==1:
                if num==0:
                    return r'$0$'
                if num==1:
                    return r'$%s$'%latex
                elif num==-1:
                    return r'$-%s$'%latex
                else:
                    return r'$%s%s$'%(num,latex)
            else:
                if num==1:
                    return r'$\frac{%s}{%s}$'%(latex,den)
                elif num==-1:
                    return r'$\frac{-%s}{%s}$'%(latex,den)
                else:
                    return r'$\frac{%s%s}{%s}$'%(num,latex,den)
        return _multiple_formatter

    class Multiple:
        def __init__(self, denominator=2, number=np.pi, latex='\pi'):
            self.denominator = denominator
            self.number = number
            self.latex = latex

        def locator(self):
            return plt.MultipleLocator(self.number / self.denominator)

        def formatter(self):
            return plt.FuncFormatter(multiple_formatter(self.denominator, self.number, self.latex))

    
    Rej_pore_array[m] = Rej_pore_inner[1:len(lit_reject)]    
    Final_rej_array[m] = np.average(R_array)
 
    
    x = Pore_array_steps
    
    
    
    fig1 = plt.figure(m)
    fig1.suptitle(comp_name_mem[m])
    ax1 = fig1.add_subplot()
    ax1.set_ylabel('Rejection (%)')
    ax1.set_xlabel('Pore radius (m)')

    
  
    
    plt.plot(Pore_array_steps, Rej_pore_array[m], label = [0, '$\pi/4$', '$\pi/2$' , '$3\pi/4$', '$\pi$', '$5\pi/4$', '$3\pi/2$', '$7\pi/4$', '$2\pi$' ])
    plt.scatter(pore_sz_plot[m], Rej_obs[m], label = 'Observed Rejection')
    plt.legend()
    ax1.get_legend().set_title("Aproach angle (Radians)")

vals = [[None for y in range(0,8)] for x in range(0,len(lit_reject))]


for i in range(0,len(lit_reject)):
    angle = list(Rej_pore_array[i][0])
    vals[i] = angle
    

max_vals = [None for x in range(0,len(lit_reject))]

for i in range(0,len(vals)):
    max_iter = vals[i]
    max_vals[i] = max(max_iter)
    
print(max_vals)

min_vals = [None for x in range(0,len(lit_reject))]

for i in range(0,len(vals)):
    min_iter = vals[i]
    min_vals[i] = min(min_iter)
    
print(min_vals)