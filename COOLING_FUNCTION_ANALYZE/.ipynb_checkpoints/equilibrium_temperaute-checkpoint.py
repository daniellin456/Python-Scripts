import numpy as np
import matplotlib.pyplot as plt

G0 = 1./1.7
Low_T_array = np.logspace(  0, 4, 1000 )
nH_array = np.logspace( -2, 11, 1000 )

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
	equilibrium_array.append( Low_T_array[index] )

plt.figure( figsize=( 20, 10 ))
plt.loglog( nH_array, equilibrium_array )
plt.plot()
plt.savefig('equilibrium_temp.pdf', dpi=600 )
plt.show()
