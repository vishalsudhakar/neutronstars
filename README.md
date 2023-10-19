# neutronstars
The library includes a 4th order Runga Kutta
method to numerical solve the Relativistic Hydrodynamical 
Equations given an Equation of State. 

Following Oppenheimer & Volkoff (1939) and Chandrasekher (1957)
the equations are re-written in terms of dimensionless quantities. I shall include 
a document illustrating the transformations.

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
