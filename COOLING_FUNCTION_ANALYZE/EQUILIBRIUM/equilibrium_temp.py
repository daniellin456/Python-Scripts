import numpy as np
import matplotlib.pyplot as plt

G0 = 1./1.7
G0 = 1.0 * G0
Low_T_array = np.logspace(  -2, 8, 10000 )
nH_array = np.logspace( -2, 8, 1000 )

def Cooling_low( nH, x, T ):

    param = G0 * np.sqrt(T) / (nH*x)
    bet = 0.74/(T**0.068)
    
    Cooling_CII   = 92.0 * 1.38E-16 * 2.0 * (( 2.8E-7 * ( T/100.0 ) ** -0.5 ) * x + 8e-10 * (T/100.0)** 0.07)* 3.5e-4 * 0.4 * np.exp(-92.0/ T)
    Cooling_O     = 1e-26 * np.sqrt(T) * ( 24.0 * np.exp( -228.0 / T) + 7.0 * np.exp( -326.0 / T ) ) * 4.5E-4
    Cooling_H     = 7.3E-19 * x * np.exp( -118400./ T ) 
    Cooling_CII_M = 6.2e4 * 1.38e-16 * 1.0 * ( 2.3e-8 * ( T / 10000.0 ) ** -0.5 * x + 1e-12 ) * np.exp( -6.2e4 / T) * 3.5e-4 * 0.4 

    Cooling_O_M = 2.3e4 * 1.38e-16 / 3.0 * ( 5.1e-9 * ( T/10000.) ** 0.57 * x + 1e-12) * np.exp(-2.3e4/T)
    Cooling_O_M = Cooling_O_M + 4.9e4 * 1.38e-16 / 3.0 * ( 2.5e-9 * (T/10000.)**(0.57) * x + 1e-12) * np.exp(-4.9e4/T)
    Cooling_O_M = Cooling_O_M + 2.6e4 * 1.38e-16 * 1.0 * ( 5.2e-9 * (T/10000.)**(0.57) * x + 1e-12) * np.exp(-2.6e4/T)
    Cooling_O_M = Cooling_O_M * 4.5e-4

    Cooling_Rec = 4.65E-30* ( T**0.94 )*( param ** bet ) * x 
    return nH ** 2 * (Cooling_CII  + Cooling_H  + Cooling_O  + Cooling_CII_M +  Cooling_CII_M + Cooling_Rec)

def Heating( nH, x, T ):
    param = G0 * np.sqrt(T) / ( nH * x ) 
    epsilon = 4.9E-2 / (1.0 + (param/1925.)**0.73) + 3.7E-2 * (T/1e4)**0.7 / (1. + (param/5e3) )
    return 1e-24 * epsilon * G0 * nH


#G0_array = [0.1*G0, 0.3*G0, 0.5*G0, 0.7*G0, G0]
#G0_factor = [0.1, 0.3, 0.5, 0.7, 1.0]
#for G0,factor in zip( G0_array, G0_factor ):
    
Total_Cooling_low = []
Total_Heating = []
for nH in nH_array:
    ne = 2.4e-3 * ( (Low_T_array/100.0) ** 0.25) / 0.5 
    x_ana = ne / nH
    x_ana = np.where( x_ana < 0.1, x_ana, 0.1 )
    x_ana = np.where( x_ana > 1.4e-4, x_ana, 1.4e-4 )
    x = x_ana
    Total_Cooling_low.append( Cooling_low( nH, x, Low_T_array ) )
    Total_Heating.append( Heating( nH, x, Low_T_array ) )
Log_Total_Cooling_low = np.log10(Total_Cooling_low)
Log_Total_Heating = np.log10(Total_Heating)
equilibrium_array = []
for i in range( 0, len( Log_Total_Cooling_low ) ):
    index = np.argsort( np.abs( np.subtract( Log_Total_Cooling_low[i], Log_Total_Heating[i] ) ) )[0]
    equilibrium_array.append(Low_T_array[index])

np.save( '/data/daniellin/PYTHON_SCRIPTS/COOLING_FUNCTION_ANALYZE/EQUILIBRIUM/nH.npy', nH_array )
np.save( '/data/daniellin/PYTHON_SCRIPTS/COOLING_FUNCTION_ANALYZE/EQUILIBRIUM/equilibrium_Temp.npy', equilibrium_array )
#print( nH_array )
#print( equilibrium_array )
#plt.figure( figsize=( 5, 4 ))
#plt.loglog( nH_array, equilibrium_array, lw=2, linestyle='-', label=str(factor)+'G0' )
#plt.tick_params(which='major',width = 1, length = 5 )
#plt.tick_params(which='minor',width = 1, length = 3 )
#plt.plot()
#plt.legend()
#plt.savefig('/data/dadiff_G0_equilibrium_temp.pdf', dpi=600 )

#plt.show()
