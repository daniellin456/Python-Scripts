import numpy as np
from  module_core import renormalize_units, normalisation


boxlen, len_pc, len_cm, vel_ms, lmax, lmin, pressure_P, temperature_K, mass_sol, mass_kg, dens_gcc,unit_J, time_Myr, mag_gauss, unit_Ek, unit_Eg0, unit_Ep, unit_Eb = renormalize_units('/drf/projets/alfven-data/phennebe/DC_1_res/')
print 'boxlen', boxlen
print 'dens_gcc', dens_gcc
print 'len_cm', len_cm
print 'vel_cms', vel_ms*100
print 'unit_Ek', unit_Ek
print 'mass_g', mass_kg*1000
#print 'Myr',time_Myr

print

pc,Ms,Myr,scale_n,scale_d,scale_t,scale_l,scale_v,scale_T2,scale_ener,scale_mag,microG,km_s,Cwnm,scale_mass,unit_col,lbox_pc = normalisation(boxlen)
print 'dens_gcc', scale_d
print 'len_cm', scale_l
print 'vel_cms', scale_v
print 'mass_g',scale_mass
print 'unit_Ek',scale_ener
#print 'Myr',Myr
#print 'Ms', Ms
