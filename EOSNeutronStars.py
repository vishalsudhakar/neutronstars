'''
This library includes a 4th order Runga Kutta
method to numerical solve the Relativistic Hydrodynamical 
Equations given an Equation of State. 

Following Oppenheimer & Volkoff (1939) and Chandrasekher (1957)
the equations are re-written in terms of dimensionless quantities.

The Equations of States (EOS) used here were derived by Wiringa et al (1988) and the fit
functions were given by Kutschera & Kotlorz (1993).

New EOS can be replaced with the ones presented here and the numerical 
calculating functions should be self consistent. It is critical, however, that the EOS are 
a function of the baryon density and the initial conditions such as initial radius, mass, and
baryon density are specified. For the EOS presented here, you can use the following initial conditions

h = (any value < 1, smaller value more precision but will run longer)
R0 = 0
t0 = 0.5 (For UV14-TN1 use t0 = 0.829)
u0 = 0

Make sure to install the following packages:

numpy
matplotlib
math

Feel free to contact me regarding any question at vsudhakar7@gatech.edu
'''

import math
import numpy as np
import matplotlib.pyplot as pt

global h_c 
global hbar_c
global m_n
global c_energydensity

# defining fundamental constants
h_c = 1240
hbar_c = h_c/(2*np.pi) # Mev fm
m_n = 939.565           # MeV/c^2

# characteristic energy density
c_energydensity = ((np.pi)**2)*((m_n)**4)/((h_c)**3)

class EOS:
    '''
    Class used to define Equations of State and its first and second derivatives
    '''

    def E(self, n, i):
        '''
        Returns the chosen EOS:
        i = 1 -> AV14+UVII
        i = 2 -> UV14-UVII
        i = 3 -> UV14-TN1

        Parameters:
        n   : the baryon density[fm^-3]
        i   : EOS number

        '''
        # input fm^-3
    
        if i == 1:
            # AV14+UVII
            func = 2.6511 + 76.744*n - 183.611*(n**2) + 459.906*(n**3) - 122.832*(n**4)   
        elif i == 2:
            # UV14-UVII
            func = 7.57891 - 1.23275*n + 227.384*(n**2) - 146.596*(n**3) + 324.823*(n**4) - 120.355*(n**5) 
        elif i == 3:
            # UV14-TN1
            func = 6.33041 - 28.1793*n + 288.397*n**2 - 65.2281*n**3
        else:
            func = 1
            
        """AV14+UVII
        func = 2.6511 + 76.744*n - 183.611*(n**2) + 459.906*(n**3) - 122.832*(n**4)
        
        UV14-UVII
        func = 7.57891 − 1.23275*n + 227.384*(n**2) − 146.596*(n**3) + 324.823*(n**4) − 120.355*(n**5) 
        
        UV14-TN1
        func = 6.33041 − 28.1793*n + 288.397*n**2 − 65.2281*n**3"""
        
        return func   # outputs MeV

    def dE_dn(self, n,i):
        '''
        Returns the chosen EOS' first derivative:
        i = 1 -> AV14+UVII
        i = 2 -> UV14-UVII
        i = 3 -> UV14-TN1

        Parameters:
        n   : the baryon density[fm^-3]
        i   : EOS number

        '''
        if i == 1:
            # AV14-UVII
            func = 76.744 - 2*(183.611)*n + 3*(459.906)*n**2 - 4*(122.832)*n**3
        elif i == 2:
            # UV14-UVII
            func = - 1.23275 + 2*227.384*n - 3*146.596*(n**2) + 4*324.823*(n**3) - 5*120.355*(n**4)
        elif i == 3:
            # UV14-TN1
            func = - 28.1793 + 2*288.397*n - 3*65.2281*n**2
        else:
            func = 1
        
        """AV14-UVII
        func = 76.744 - 2*(183.611)*n + 3*(459.906)*n**2 - 4*(122.832)*n**3
        
        UV14-UVII
        unc = − 1.23275 + 2*227.384*n − 3*146.596*(n**2) + 4*324.823*(n**3) − 5*120.355*(n**4) 
        
        UV14-TN1
        func = − 28.1793 + 2*288.397*n − 3*65.2281*n**2""" 
        
        return func   # outputs MeV * fm^-3

    def d2E_dn2(self, n,i):
        '''
        Returns the chosen EOS' second derivative:
        i = 1 -> AV14+UVII
        i = 2 -> UV14-UVII
        i = 3 -> UV14-TN1

        Parameters:
        n   : the baryon density[fm^-3]
        i   : EOS number

        '''
        
        if i == 1:
            # AV14-UVII
            func = -2*(183.611) + 3*2*(459.906)*n - 4*3*(122.832)*n**2
        elif i == 2:
            # UV14-UVII
            func = 2*227.384 - 2*3*146.596*n + 3*4*324.823*(n**2) - 4*5*120.355*(n**3)
        elif i == 3:
            # UV14-TN1
            func = 2*288.397 - 2*3*65.2281*n
        else:
            func = 1
            
        """AV14-UVII
        func = -2*(183.611) + 3*2*(459.906)*n - 4*3*(122.832)*n**2
        
        UV14-UVII
        func = 2*227.384 − 2*3*146.596*n + 3*4*324.823*(n**2) − 4*5*120.355*(n**3) 
        
        UV14-TN1
        func = 2*288.397 − 2*3*65.2281*n"""
        
        return func  # outputs MeV * fm^-6
    

class HydrodynamicalEquations():
    '''
    The class encoperates the Relativistic Hydrodynamics Equations
    and variable transformations.
    '''

    def __init__(self):
    
        self.eos = EOS()        # EOS object

    def n(self, t):
        '''
        The baryon density as a function of unitless parameter t.
        '''
        term1 = (np.sinh(t/4))**3
        term2 = (3*(np.pi**2)*((hbar_c/(m_n))**3))
        nvalue = term1/term2

        return nvalue # fm**-3
    
    def pressure(self, n, i):
        '''
        outputs the pressure given the baryon density for an EOS

        n   : the baryon density[fm^-3]
        i   : EOS number
        '''
        pvalue = (n**2)*(self.eos.dE_dn(n,i))
        
        return pvalue/c_energydensity 

    def energydensity(self, n, i):
        '''
        outputs the energy density given the baryon density for an EOS

        n   : the baryon density[fm^-3]
        i   : EOS number
        '''
        evalue = n*(self.eos.E(n,i) + m_n)

        return evalue/c_energydensity  
    
    def dt_dP(self, t, i):
        '''
        outputs the dimensionless baryon density parameter as a change in presseure.
        This is required to dt/dR

        t   : the dimensionless baryon density parameter
        i   : EOS number
        '''
        nvalue = self.n(t)
        
        term1 = 2*(nvalue)*self.eos.dE_dn(nvalue,i) + (nvalue**2)*(self.eos.d2E_dn2(nvalue,i))
        term2 = 3*((np.sinh(t/4))**2)*(np.cosh(t/4))/(12*((np.pi)**2)*(hbar_c/(m_n))**3)
        
        value = (term1*term2)/c_energydensity
        
        if value != 0:
            return value**-1
        else:
            print("dp_dt = 0")
            return 0 
        
    
    def dt_dR(self, R, t, u, i):
        '''
        returns the dt/dR for a particular dimensionless radius R

        Parameters:
        R   : dimensionless radius
        t   : dimensionless baryon density
        u   : dimensionless mass
        i   : EOS number
        '''
        term1 = -4*np.pi*R
        term2 = self.dt_dP(t, i)
        term3 = (self.pressure(self.n(t),i) + self.energydensity(self.n(t),i))/(1 - (2*u/R))
        term4 = self.pressure(self.n(t),i) + (u/(4*np.pi*R**2))
        
        value = term1*term2*term3*term4
        return value  
    
    def dt_dR1(self, R, t, i):
        '''
        returns the dt/dR when we reach the core of the neutron star. This is to
        ensure we don't hit singularities.

        Parameters:
        R   : dimensionless radius
        t   : dimensionless baryon density
        i   : EOS number
        '''
        term1 = -4*np.pi*R
        term2 = self.dt_dP(t,i)
        term3 = (self.pressure(self.n(t),i) + self.energydensity(self.n(t),i))
        term4 = self.pressure(self.n(t),i)
        
        value = term1*term2*term3*term4

        return value 

    
    def du_dR(self, R, t, i):
        '''
        returns the du/dR when we reach the core of the neutron star

        Parameters:
        R   : dimensionless radius
        t   : dimensionless baryon density
        i   : EOS number
        '''
        derivative = 4*(np.pi)*self.energydensity(self.n(t),i)*R**2

        return derivative 
    
    def runga(self, h, R0, t0, u0, i):
        '''
        performs the 4th order Runga Kutta, increasing the radius by h
        until the baryon density and pressure greater than and equal 
        to 0. Returns the mass and radius of the neutron star for a particular
        baryon density.

        Parameters:
        h   : the step in radium delta R
        R0  : initial radius starting point
        t0  : initial baryon density starting point
        u0  : initial mass starting point
        i   : EOS number
        '''
        un = u0
        tn = t0
        R = R0
        
        Ri = R
        ui = un
        
        nv = self.n(tn)
        pv = self.pressure(nv,i)
        while(nv >= 0 and pv >= 0):
            ui = un
            ti = tn
            Ri = R

            k1_u = h*self.du_dR(R,ti,i)
            if(R == 0):
                k1_t = h*self.dt_dR1(R,ti,i)
            else:
                k1_t = h*self.dt_dR(R,ti,ui,i)

            # check n and p
            nv = self.n(ti + k1_t/2)
            pv = self.pressure(nv, i)

            if(nv < 0 or pv <0):
                break


            k2_u = h*self.du_dR(R + h/2, ti + k1_t/2, i)
            k2_t = h*self.dt_dR(R + h/2, ti + k1_t/2, ui + k1_u/2, i)

            # check n and p
            nv = self.n(ti + k2_t/2)
            pv = self.pressure(nv,i)
            
            if(nv < 0 or pv < 0):
                break

            k3_u = h*self.du_dR(R + h/2, ti + k2_t/2, i)
            k3_t = h*self.dt_dR(R + h/2, ti + k2_t/2, ui + k2_u/2, i)

            # check n and p
            nv = self.n(ti + k3_t)
            pv = self.pressure(nv,i)

            if(nv < 0 or pv < 0):
                break

            k4_u = h*self.du_dR(R + h,ti + k3_t,i)
            k4_t = h*self.dt_dR(R + h,ti + k3_t, ui + k3_u,i)

            # check n and p
            tn = ti + k1_t/6 + k2_t/3 + k3_t/3 + k4_t/6
            un = ui + k1_u/6 + k2_u/3 + k3_u/3 + k4_u/6
            R = Ri + h        

            nv = self.n(tn)
            pv = self.pressure(nv,i)
        
        return (self.n(t0), (Ri*13.69), ui*9.29)
    
class run(EOS, HydrodynamicalEquations):
    '''
    contains the run functions and plot functions that get the mass and radius for a range
    of baryon densities.
    '''
    def __init__(self):
        super().__init__()
        
        self.cd = np.array([])
        self.r = np.array([])
        self.m = np.array([])

        self.complete_data = ()

        pt.rcParams['font.weight'] = 'normal'
        pt.rcParams['axes.labelweight'] = 'normal'
        pt.rcParams['axes.linewidth'] = 1.5
        pt.rcParams['lines.linewidth'] = 1.5


    def loop(self, h, R0, t0, u0, i):
        '''
        performs the numerical calculations for a range of baryon 
        density and returns three arrays representing 
        the (baryon density, radius, mass)

        Parameters:

        '''
        # initial conditions 
        R = R0
        t = t0
        u = u0
        
        # data set
        cd = []
        r = []
        mass = []
        
        # stop condition
        dP_dt = (self.dt_dP(t,i))**-1
        
        while (dP_dt >= 0):
            
            each_data = self.runga(h,R,t,u,i)
            cd.append(each_data[0])
            r.append(each_data[1])
            mass.append(each_data[2])
            t = t + 0.01
            print("Calculating for baryon density t = " + str(t))
            dP_dt = (self.dt_dP(t,i))**-1
            print("")
            
        return (cd, r, mass)

    def calculate(self, h, R0, t0, u0, i):
        
        data = self.loop(h, R0, t0, u0, i)

        self.cd = np.array(data[0])
        self.r = np.array(data[1])
        self.m = np.array(data[2])

        self.complete_data = (self.cd, self.r, self.m)

        
    def plot_mass_baryon_density(self, label, color = (173/255, 52/255, 62/255), savefig=False):
        '''
        creates a plot of the mass vs the baryon density

        Parameters:
        label (str)         : the label of the curve
        color (str)         : the color of the curve
        savefig (boolean)   : boolean value to save the plot
        '''

        pt.plot(self.cd, self.m, label = label, color=color)
        pt.ylabel("Mass ($M_{\odot}$)", fontsize=13)
        pt.xlabel("n$_{0}$ (fm$^{-3}$)", fontsize=13)
        pt.xlim([0,2.5])
        pt.ylim([0,3])
        pt.legend(loc=4, frameon=False, fontsize=13)
        pt.minorticks_on()
        pt.tick_params(which='minor', direction='in', width=2)
        pt.tick_params(which='major', direction='in', width=2)
        if savefig:
            pt.savefig("mass_n0.pdf")
        pt.show()

    def plot_radius_baryon_density(self, label='', color = (173/255, 52/255, 62/255), savefig=False):
        '''
        creates a plot of the radius vs baryon density

        Parameters:
        label (str)         : the label of the curve
        color (str)         : the color of the curve
        savefig (boolean)   : boolean value to save the plot
        '''

        pt.plot(self.cd, self.r, label = label, color=color)
        pt.ylabel("Radius (km)", fontsize=13)
        pt.xlabel("n$_{0}$ (fm$^{-3}$)", fontweight='bold', fontsize=13)
        pt.xlim([0,2.5])
        pt.ylim([0,15])
        pt.legend(loc=4, frameon=False, fontsize=13)
        pt.minorticks_on()
        pt.tick_params(which='minor', direction='in', width=2)
        pt.tick_params(which='major', direction='in', width=2)
        if savefig:
            pt.savefig("radius_n0.pdf")
        pt.show()

    def plot_mass_radius(self, label='', color = (173/255, 52/255, 62/255), savefig=False):
        '''
        creates a plot of the mass vs radius

        Parameters:
        label (str)         : the label of the curve
        color (str)         : the color of the curve
        savefig (boolean)   : boolean value to save the plot
        '''

        pt.plot(self.r, self.m, label = label, color=color)
        pt.ylabel("Mass (M$_{\odot}$)", fontsize=13)
        pt.xlabel("Radius (km)", fontsize=13)
        pt.xlim([6,20])
        pt.ylim([0,3])
        point_x = [13.02,13.7]
        point_y = [1.44, 2.08]
        error_x = [[1.06,1.24],[1.5,2.6]]
        error_y = [[0.14,0.07],[0.15,0.07]]
        pt.annotate("PSR J0030+0451",(13.02,1.44),textcoords="offset points",xytext=(10,10), fontsize=11)
        pt.annotate("PSR J0740+6620",(13.7,2.08),textcoords="offset points",xytext=(10,10), fontsize=11)
        pt.errorbar(point_x,point_y, xerr = error_x, yerr = error_y, fmt = "o",color = "black")
        pt.legend(loc=4, frameon=False, fontsize=11)
        pt.minorticks_on()
        pt.tick_params(which='minor', direction='in', width=2)
        pt.tick_params(which='major', direction='in', width=2)
        
        if savefig:
            pt.savefig("mass_radius.pdf")
            
        pt.show()

