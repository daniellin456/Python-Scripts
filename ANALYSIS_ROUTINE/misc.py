import argparse
import numpy as np
import matplotlib
#matplotlib.use('agg')
import matplotlib.pyplot as plt
from pymses.filters import CellsToPoints
from pymses import RamsesOutput
from pymses.sources.ramses.filename_utils import search_valid_outputs as search_ro
import sys

def read_amr(readdir,ioutput=None):
	
	outputs = search_ro(readdir)
	print outputs

	# select samples 
	#nout = len(outputs)
	#isample = nout/3
	#if(isample==0):isample=1
	#outputs=outputs[-1:-1-3*isample:-isample]
	#nout = len(outputs)
	#for i in range(nout)[:1]:
	#    ioutput=outputs[i]


	if ioutput is None: ioutput=outputs[-1]
	print ioutput
	return RamsesOutput(readdir, ioutput)

def get_data(amr, dens, vol_cc,  weight="volume"):
	assert weight in ["mass", "volume", "count"]

	cell_src = CellsToPoints(amr)#,smallest_cell_level=8)
	cells = cell_src.flatten()

	dx = cells.get_sizes()
	vol = dx**3. * vol_cc


	data = cells["rho"]*dens
	weights = None
	if weight == "mass":
		weights = cells["rho"] * vol * dens_gcm
	elif weight == "volume":
		weights = vol
	else:
		weights = np.ones_like(vol)

	return data, weights

def get_data_cells(cells, dens, vol_cc,  weight="volume"):
    assert weight in ["mass", "volume", "count"]

    dx = cells.get_sizes()
    vol = dx**3. * vol_cc

    data = cells["rho"]*dens
    weights = None
    if weight == "mass":
        weights = cells["rho"] * vol * dens_gcm
    elif weight == "volume":
        weights = vol
    else:
        weights = np.ones_like(vol)

    return data, weights

def make_bins(data, logx, nbins):

    bins = None

    a = data.min() / 1.01
    b = data.max() * 1.01
    if logx:
        bins = a * (b / a) ** (np.arange(nbins + 1, dtype=np.float) / float(nbins))
    else:
        bins = a + (b - a) * (np.arange(nbins + 1, dtype=np.float) / float(nbins))
    return bins


def get_hist(data, logx, logy, nbins, normalize, weight=np.array([]), plot=False, color='b',lw=1,label=None):

	# data: numpy.ndarray

	bins = make_bins(data, logx, nbins)

	if weight.any():
		y, plot_bins = np.histogram(data, bins, weights=weight, density=normalize)
	else:
		y, plot_bins = np.histogram(data, bins, density=normalize)

	y=y.astype(float)


	y = y / np.sum(y) / np.diff(bins)

	if logx: 
		x = np.sqrt(bins[:-1]*bins[1:])
		y=y*x
	else: x=(bins[:-1]+bins[1:])/2.


	if plot == True:
		plots = {(False, False): plt.plot,
				(False,  True): plt.semilogy,
				( True, False): plt.semilogx,
				( True,  True): plt.loglog}

		plots[logx, logy](x, y, color,lw=lw,label=label)

	return x,y


def get_leaf_hist(data, data_l,logx, logy, nbins, normalize, weight=np.array([]), plot=False, color='b',lw=1,label=None):

	# data: numpy.ndarray

	bins = make_bins(data, logx, nbins)

	if weight.any():
		y, plot_bins = np.histogram(data, bins, weights=weight, density=normalize)
		y_l, plot_bins = np.histogram(data_l, bins, weights=weight, density=normalize)
	else:
		y, plot_bins = np.histogram(data, bins, density=normalize)
		y_l, plot_bins = np.histogram(data_l, bins, density=normalize)

	y=y.astype(float)

	y_l = y_l / np.sum(y) / np.diff(bins)

	if logx: 
		x = np.sqrt(bins[:-1]*bins[1:])
		y_l=y_l*x
	else: x=(bins[:-1]+bins[1:])/2.


	if plot == True:
		plots = {(False, False): plt.plot,
				(False,  True): plt.semilogy,
				( True, False): plt.semilogx,
				( True,  True): plt.loglog}

		plots[logx, logy](x, y_l, color,lw=lw,label=label)

	return x,y_l


def get_cdf_hist(data, logx, logy, nbins, normalize, weight=np.array([]), plot=False, color='b',lw=1,label=None):

	# data: numpy.ndarray

	bins = make_bins(data, logx, nbins)

	if weight.any():
		y, plot_bins = np.histogram(data, bins, weights=weight, density=normalize)
	else:
		y, plot_bins = np.histogram(data, bins, density=normalize)

	y=y.astype(float)
	y = y / np.sum(y) / np.diff(bins)

	if logx: 
		x = np.sqrt(bins[:-1]*bins[1:])
	else: x=(bins[:-1]+bins[1:])/2.

	cdf = np.empty([len(y)])

	for i in range(len(cdf)):
		cdf[i] = (y*np.diff(bins))[:i+1].sum()

	if plot == True:
		plots = {(False, False): plt.plot,
				(False,  True): plt.semilogy,
				( True, False): plt.semilogx,
				( True,  True): plt.loglog}

		plots[logx, logy](x, cdf, color,lw=lw,label=label)

	return x,cdf

def check_zero(data):

    for i in range(len(data)):
        if i==0 :
            if data[i] == 0: data[i]=data[i+1]
        elif i==(len(data)-1) :
            if data[i] == 0: data[i]=data[i-1]
        else:
            if data[i] == 0: data[i]=(data[i-1]+data[i+1])/2

    return data


__all__ = ["read_amr","get_data","get_data_cells","get_hist","get_leaf_hist","get_cdf_hist","check_zero"]
