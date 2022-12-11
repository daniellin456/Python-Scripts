import numpy as np
import matplotlib
from scipy import special

def MLP (m,mu,std,alpha):

    x=(1/2**0.5)*(alpha*std - (np.log(m)-mu)/std)

    return (alpha/2)*np.exp(alpha*mu+alpha**2*std**2/2)*m**(-1-alpha)*special.erfc(x)

def mMLP (m,mu,std,alpha):

    x=(1/2**0.5)*(alpha*std - (np.log(m)-mu)/std)

    return m*(alpha/2)*np.exp(alpha*mu+alpha**2*std**2/2)*m**(-1-alpha)*special.erfc(x)

def logMLP (m,mu,std,alpha):

    x=(1/2**0.5)*(alpha*std - (np.log(m)-mu)/std)

    y=(alpha/2)*np.exp(alpha*mu+alpha**2*std**2/2)*m**(-1-alpha)*special.erfc(x)
    return np.log10(y)

def logmMLP (m,mu,std,alpha):

	x=(1/2**0.5)*(alpha*std - (np.log(m)-mu)/std)

	y =  m*(alpha/2)*np.exp(alpha*mu+alpha**2*std**2/2)*m**(-1-alpha)*special.erfc(x)
	return np.log10(y)

def MLP_cdf (m,mu,std,alpha):

	x=(1/2**0.5)*(alpha*std - (np.log(m)-mu)/std)

	x2=-(np.log(m)-mu)/(2**0.5*std)

	return (1./2.)*special.erfc(x2)-(1./2.)*np.exp(alpha*mu+alpha**2*std**2/2)*m**(-alpha)*special.erfc(x)

def lognormal(x,std):

	mu =-std**2/2

	return (1/std/(2*np.pi)**0.5)*np.exp(-(np.log(x)-mu)**2/(2*std**2))

def logoflognormal(x,std,mu):

	y= (1/x/std/(2*np.pi)**0.5)*np.exp(-(np.log(x)-mu)**2/(2*std**2))
	return np.log10(y)


__all__ = ["MLP","mMLP","logMLP","logmMLP","MLP_cdf","lognormal","logoflognormal"]
