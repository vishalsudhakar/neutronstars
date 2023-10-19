'''
This file illustrates how to use the numerical
calculation library, eos_neutron_stars.py, to get
mass and radius of neutron star as a function of the 
baryon density for AV14-UVII EOS.

Make sure to install the following packages:

numpy
matplotlib
math
'''
import EOSNeutronStars as eos

main = eos.run()

# specify initial condition (notes these are dimensionless parameters)
h = 0.0007            # radius step
R0 = 0               # starting radius
t0 = 0.5             # starting baryon density
u0 = 0               # starting mass
i  = 1               # choosing the EOS (AV14-UVII)

main.calculate(h, R0, t0, u0, i)
main.plot_mass_baryon_density(label="AV14-UVII")
main.plot_radius_baryon_density(label="AV14-UVII")
main.plot_mass_radius(label="AV14-UVII")


