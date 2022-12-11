import numpy as np
import matplotlib.pyplot as plt

#G0 = 1./1.7
G0 = 0.0
Low_T_array = np.logspace(  1, 4, 100 )
#High_T_array = np.logspace(  4.1, 7, 100 )
nH_array = np.logspace( -2, 6, 100 )
Color = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']


def Cooling_high( nH, x, T ):
    logT = np.log10(T)
    if (logT < 4.0):
        Cooling = 0.1343 * logT ** 3 - 1.3906 * logT ** 2 + 5.1554 * logT - 31.967
    elif (logT < 4.25):
        Cooling = 12.64 * logT - 75.56
    elif (logT < 4.35):
        Cooling = -0.3 * logT - 20.565
    elif (logT < 4.9):
        Cooling = 1.745 * logT - 29.463
    elif (logT < 5.4):
        Cooling = -20.9125
    elif (logT < 5.9):
        Cooling = -1.795 * logT - 11.219
    elif (logT < 6.2):
        Cooling = - 21.8095
    elif (logT < 6.7):
        Cooling = -1.261 * logT -13.991
    else:
        Cooling = -22.44

    return -1.0*10.0**(Cooling)

def Cooling_low( nH, x, T ):

    param = G0 * np.sqrt(T) / (nH*x)
    bet = 0.74/(T**0.068)
    
    Cooling_CII   = 92.0 * 1.38E-16 * 2.0 * (( 2.8E-7 * ( T/100.0 ) ** -0.5 ) * x + 8e-10 * (T/100.0)** 0.07)* 3.5e-4 * 0.4 * np.exp(-92.0/ T)
    Cooling_O     = 1e-26 * np.sqrt(T) * ( 24.0 * np.exp( -228.0 / T) + 7.0 * np.exp( -326.0 / T ) ) * 4.5E-4
    Cooling_H     = 7.3E-19 * x * np.exp( -118400./ T )
    Cooling_CII_M = 6.2e4 * 1.38e-16 * 1.0 * ( 2.3e-8 * ( T / 10000.0 ) ** -0.5 * x + 1e-12 ) * np.exp( -6.2e4 / T) * 3.5e-4 * 0.4

    #if ( T < 1e4 ):
    Cooling_O_M = 2.3e4 * 1.38e-16 / 3.0 * ( 5.1e-9 * ( T/10000.) ** 0.57 * x + 1e-12) * np.exp(-2.3e4/T)
    Cooling_O_M = Cooling_O_M + 4.9e4 * 1.38e-16 / 3.0 * ( 2.5e-9 * (T/10000.)**(0.57) * x + 1e-12) * np.exp(-4.9e4/T)
    Cooling_O_M = Cooling_O_M + 2.6e4 * 1.38e-16 * 1.0 * ( 5.2e-9 * (T/10000.)**(0.57) * x + 1e-12) * np.exp(-2.6e4/T)
    #else:
    #    Cooling_O_M = 2.3e4 * 1.38e-16 / 3.0 * ( 5.1e-9 * (T/10000.) ** 0.17 * x + 1e-12) * np.exp(-2.3e4/T)
    #    Cooling_O_M = Cooling_O_M + 4.9e4 * 1.38e-16 / 3.0  * ( 2.5e-9 * (T/10000.) ** 0.13 * x + 1e-12) * np.exp(-4.9e4/T)
    #    Cooling_O_M = Cooling_O_M + 2.6e4 * 1.38e-16 * 1.0  * ( 5.2e-9 * (T/10000.) ** 0.15 * x + 1e-12) * np.exp(-2.6e4/T)
    
    Cooling_O_M = Cooling_O_M * 4.5e-4
    Cooling_Rec = 4.65E-30* ( T**0.94 )*( param ** bet ) * x
    return nH ** 2 * (Cooling_CII  + Cooling_H  + Cooling_O  + Cooling_CII_M +  Cooling_CII_M + Cooling_Rec)

def Heating( nH, x, T ):
    param = G0 * np.sqrt(T) / ( nH * x )
    epsilon = 4.9E-2 / (1.0 + (param/1925.)**0.73) + 3.7E-2 * (T/1e4)**0.7 / (1. + (param/5e3) )
    return 1e-24 * epsilon * G0 * nH

Total_Cooling_high = []
Total_Cooling_low = []
Total_Heating = []

for nH in nH_array:
    #for T in T_array:
        ne = 2.4e-3 * ( (Low_T_array/100.0) ** 0.25) / 0.5
        x_ana = ne / nH
        x_ana = np.where( x_ana < 0.1, x_ana, 0.1 )
        x_ana = np.where( x_ana > 1.4e-4, x_ana, 1.4e-4 )
        x = x_ana
        #if ( T > 10035 ):
        #    Total_Cooling.append( Cooling_high( nH, x, T ) )
        #    Total_Heating.append(0.0)
        #else:
        Total_Cooling_low.append( Cooling_low( nH, x, Low_T_array ) )
        Total_Heating.append( Heating( nH, x, Low_T_array ) )
        
#for nH in nH_array:
#    for T in High_T_array:
#        Total_Cooling_high.append( Cooling_high( nH, x, T ) )


plt.figure( figsize=(20, 10) )
for i in range( 0, 100, 10 ):
    #plt.loglog( High_T_array, np.fabs(Total_Cooling_high[ i*100 : (i+1)*100]), Color[int(i/10)] )
    plt.loglog( Low_T_array, Total_Cooling_low[i], Color[int(i/10)], label='nH=' + str(nH_array[i]) )
    plt.loglog( Low_T_array, Total_Heating[i], Color[int(i/10)], linestyle='--', linewidth=2 )
plt.legend()
plt.savefig('ramses_cooling_function.pdf', dpi=600)
