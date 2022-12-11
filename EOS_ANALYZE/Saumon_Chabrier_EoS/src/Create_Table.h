#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

double EoS_DensEint_Temp( const double Dens, const double Eint );
double EoS_DensEint_Pres( const double Dens, const double Eint );
void Construct_DensTemp2Eint_Table();
void Construct_DensPres2Eint_Table();
