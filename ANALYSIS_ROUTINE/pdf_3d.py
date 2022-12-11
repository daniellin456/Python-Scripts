#! /usr/bin/env python

import argparse
import numpy as np
import matplotlib
#matplotlib.use('agg')
import matplotlib.pyplot as plt
from pymses.utils import constants as C
from scipy.optimize import curve_fit
from pymses.utils.regions import Sphere, SphericalShell
import pymses
from pymses.filters import CellsToPoints
from pymses import RamsesOutput
from pymses.sources.ramses.filename_utils import search_valid_outputs as search_ro
from pymses.filters import RegionFilter#, RegionDifferentiator
from pymses.sources.ramses.amr	import read_ramses_amr_file
#from pymses.utils.regions import Ellipsoid
from misc import *
from func import *
import os,sys


if __name__ == '__main__':


	parser = argparse.ArgumentParser(description="Calculate profiles")
	parser.add_argument("out_dir", type=str, help="RAMSES output repository")
	parser.add_argument("iout",type=int, help="ioutput")
	parser.add_argument("max_read_level",type=int, help="max read level")
	args = parser.parse_args()

	readdir = '/data/ynlee/'+args.out_dir
	ro = read_amr(readdir,args.iout)


	colors = ['b','g','r','c','m','y','k']
	vol_cc = 1.#(ro.info["unit_length"].express(C.cm))**3

	dens_gcm = ro.info["unit_density"].express(C.g/C.cm**3) 
	dens_Msun = ro.info["unit_density"].express(C.Msun/C.cm**3)
	dens_Hcc = 1.#ro.info["unit_density"].express(C.H_cc)
	len_cm = ro.info["unit_length"].express(C.cm)/ro.info["boxlen"]


	#### read amr data ####
	amr = ro.amr_source(["rho"])
	#ro.amr_fields()
	amr.set_read_levelmax(args.max_read_level)

	#cells = CellsToPoints(amr).flatten()
	#print np.min(cells["rho"])

	'''#### filter data with a region ####
	center = [0.5, 0.5, 0.5]
	#radius = 0.25
	#region = Sphere(center, radius)
	r_in = 0.4
	r_out = 0.5
	region = SphericalShell(center, r_in, r_out)
	amr = RegionFilter(region,amr)'''

	#### plot histogram ####
	data, weights = get_data(amr, dens_Hcc, vol_cc, weight="volume")
	x, y = get_hist(data, logx=True, logy=True, nbins=100, normalize=False, weight=weights,plot=True,label='pdf')


#	#### fit ####
#	fit_min=1.e0
#	fit_max=1.e15
#	ind = (x<fit_max) * (x>fit_min)
#	y2 = np.log10(check_zero(y))
#
#	ini_value = np.array([4,2,1])
#
#	popt, pcov = curve_fit(logmMLP,x[ind],y2[ind],ini_value)
#
#	mu = popt[0]
#	std = popt[1]
#	alpha = popt[2]
#	slope = alpha+1
#	print 'slope',slope
#	print 'mu:',mu,'std:',std,'alpha:',alpha
#
#	plt.loglog(x[ind],mMLP(x[ind] ,popt[0],popt[1],popt[2]))
#	plt.text (1e2,1e-8,"slope:%f"%slope, fontsize=12)
#
#	# calculate transition  point
#
#	th_x = (x[ind][np.argmax(y[ind])])*10** std
#	print 'transiton point:',th_x
#	diff = abs(x[0]-th_x)
#	for i  in range(len(x)):
#		if diff > abs(x[i]-th_x):
#			th_i = i
#			diff = abs(x[i]-th_x)
#
#
#   
#	#th_y = y[ind][np.argmax(y[ind])]+ std
#	plt.loglog(x[th_i], y[th_i], marker='o', markersize=3)
#	plt.text (x[th_i],y[th_i],"threshold")
#
#	x_mass, y_mass = get_hist(data, logx=True, logy=True, nbins=100, normalize=False, plot=False, color='y--',label='2D pdf')
#	mass_tot = np.trapz(y_mass,x_mass) 
#	mass_th = np.trapz(y_mass[th_i:],x_mass[th_i:]) 
#	ratio = mass_th/mass_tot
#
#	print 'tot mass:',mass_tot
#	print 'dense mass:',mass_th
#	print 'ratio',ratio
#	#plt.text (1e8,1e-6,"total mass:%5e"%mass_tot)
#	#plt.text (1e8,1e-7,"dense mass:%5e"%mass_th)
#	plt.text (1e2,1e-9,"ratio:%5.4f"%ratio,fontsize=12)

	#plt.xticks(10.**np.arange(0,11,2),fontsize=16)
	#plt.yticks(10.**np.arange(-12,1,4),fontsize=16)
	#plt.xlim([1.e-1,1.e10])
	#plt.ylim(top=1.e-4)
	plt.xlabel(r'$\rho$ (H/cc)',fontsize=14)
	plt.ylabel(r'$d(V/V_\mathrm{tot})/d\mathrm{log}\rho$',fontsize=14)
	plt.legend()
	plt.suptitle("Read level:%3d"%args.max_read_level)



	
	result_fig = readdir+"/pdf_compact_%i_%03d"%(args.iout,args.max_read_level)+".eps"
	plt.savefig(result_fig)
	plt.show()

