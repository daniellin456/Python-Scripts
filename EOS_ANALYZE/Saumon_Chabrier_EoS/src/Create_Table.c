#include "Create_Table.h"

#define LOG10( a ) log10(a)
#define FABS( a )  fabs(a)
#define POW( a, b )   pow(a, b)
#define MAX( a, b )     (  ( (a) > (b) ) ? (a) : (b)  )
#define MIN( a, b )     (  ( (a) < (b) ) ? (a) : (b)  )
#define SQR( a ) ( ( a ) * ( a ) )
#define FLOOR(a) floor(a)

extern int nRho;
extern int nEnergy;
extern int nTemp;
extern int nPres;

extern double UNIT_D;
extern double UNIT_V;
extern double Const_kB;
extern double Const_amu;
extern double GAMMA;
extern double MOLECULAR_WEIGHT;

extern double RhoMin;
extern double RhoMax;
extern double Emin;
extern double Emax;
extern double TempMin;
extern double TempMax;
extern double PresMin;
extern double PresMax;
extern double dRho;
extern double dEner;
extern double dTemp;
extern double dPres;

extern double *Dens_Table;
extern double *Ener_Table;
extern double *DensEint2Pres_Table;
extern double *DensEint2Temp_Table;
extern double *DensEint2Cs_Table;
extern double *DensPres2Eint_Table;
extern double *DensTemp2Eint_Table;

double EoS_DensEint2Temp( const double Dens, const double Eint )
{
	if ( Dens <= 0.0 || Eint <= 0.0 ) return 0.0;
	const double Dens_CGS = Dens * UNIT_D;
	const double Eint_CGS = Eint * UNIT_D * SQR( UNIT_V );
	const double lr = (LOG10(Dens_CGS) - RhoMin) / dRho;
	const double le = (LOG10(Eint_CGS) - Emin - LOG10(Dens_CGS) ) / dEner;
	if ( lr < 0 || lr > nRho - 1 ) fprintf(stderr, "Interpolation Out of Bound, Rho = %20.14e\n", Dens_CGS );
	if ( le < 0 || le > nEnergy - 1 ) fprintf(stderr, "Interpolation Out of Bound, Eint = %20.14e\n", Eint_CGS );

	const int ir = FLOOR(lr);
	const int ie = FLOOR(le);
	const double x = lr - ir;
	const double y = le - ie;

	double Temp = 0.0;
	Temp += ( 1.0 - x ) * ( 1.0 - y ) * DensEint2Temp_Table[ ie * nRho + ir ];
	Temp +=       x     * ( 1.0 - y ) * DensEint2Temp_Table[ ie * nRho + ir + 1 ];
	Temp += ( 1.0 - x ) *         y   * DensEint2Temp_Table[ (ie + 1)* nRho + ir ];
	Temp +=       x     *         y   * DensEint2Temp_Table[ (ie + 1)* nRho + ir + 1];

	if ( DensEint2Temp_Table[ ie * nRho + ir ] *
	     DensEint2Temp_Table[ ie * nRho + ir + 1 ] *
	     DensEint2Temp_Table[ (ie + 1)* nRho + ir ] *
	     DensEint2Temp_Table[ (ie + 1)* nRho + ir + 1] == 0.0 )
	{
	   Temp = 0.0;
	}

	return Temp;
}

double EoS_DensEint2Pres( const double Dens, const double Eint )
{
   if ( Eint <= 0.0 || Dens <= 0.0 ) return 0.0;

   const double Dens_CGS = Dens * UNIT_D;
   const double Eint_CGS = Eint * UNIT_D * SQR( UNIT_V );
   const double lr = (LOG10(Dens_CGS) - RhoMin) / dRho;
   const double le = (LOG10(Eint_CGS) - Emin - LOG10(Dens_CGS)) / dEner;
   if ( lr < 0 || lr > nRho - 1 ) return 0.0; //fprintf(stderr, "Interpolation Out of Bound, Rho = %20.14e\n", Dens_CGS );
   if ( le < 0 || le > nEnergy - 1 ) return 0.0; //fprintf(stderr, "Interpolation Out of Bound, Rho = %20.14e, Eint = %20.14e, MinEner = %20.14e, dEner = %20.14e\n", Dens_CGS, Eint_CGS, Min_Ener, d_Ener );

   const int ir = FLOOR(lr);
   const int ie = FLOOR(le);
   const double x = lr - ir;
   const double y = le - ie;

   double Pres_CGS = 0.0;
   double Pres_Code = 0.0;
   Pres_CGS += ( 1 - x ) * ( 1 - y ) * DensEint2Pres_Table[ ie * nRho + ir ];
   Pres_CGS +=       x   * ( 1 - y ) * DensEint2Pres_Table[ ie * nRho + ir + 1 ];
   Pres_CGS += ( 1 - x ) *       y   * DensEint2Pres_Table[ (ie + 1)* nRho + ir ];
   Pres_CGS +=       x   *       y   * DensEint2Pres_Table[ (ie + 1)* nRho + ir + 1];
   Pres_Code = Pres_CGS / (UNIT_D * SQR(UNIT_V) );

   if ( DensEint2Pres_Table[ ie * nRho + ir ] *
        DensEint2Pres_Table[ ie * nRho + ir + 1 ] *
        DensEint2Pres_Table[ (ie + 1)* nRho + ir ] *
        DensEint2Pres_Table[ (ie + 1)* nRho + ir + 1] == 0.0 )
   {
      Pres_Code = 0.0;
   }

	return Pres_Code;
}

void Construct_DensTemp2Eint_Table()
{
    fprintf(stdout, "%s ... \n", __FUNCTION__);

    DensTemp2Eint_Table = (double *)malloc(sizeof(double) * nRho * nTemp);
    memset(DensTemp2Eint_Table, 0, sizeof(double) * nRho * nTemp);

    for (int ir = 1; ir < nRho - 1; ir++)
    {
        for (int it = 0; it < nTemp; it++)
        {
            double Rho0 = POW(10.0, Dens_Table[1 * nRho + ir]);
            double Temp0 = POW(10.0, (LOG10(TempMin) + it * dTemp));
            double Eint_Old = Rho0 * Const_kB * Temp0 / (MOLECULAR_WEIGHT * Const_amu * (GAMMA - 1.0));
            if (it > 0) Eint_Old = MAX(Eint_Old, DensTemp2Eint_Table[ ( it - 1 )* nRho + ir]);

            for (int ii = 0; ii < 1000; ii++)
            {
                double Temp1 = EoS_DensEint2Temp(Rho0/UNIT_D, Eint_Old/(UNIT_D*SQR(UNIT_V)));
                if (Temp1 == 0.0)
                {
                    Eint_Old = 0.0;
                    break;
                }

                double Temp2 = EoS_DensEint2Temp(Rho0/UNIT_D, Eint_Old*1.001/(UNIT_D*SQR(UNIT_V)));
                if (Temp2 == 0.0)
                {
                    Eint_Old = 0.0;
                    break;
                }

                double Eint_New;
                if (FABS(Temp2 - Temp1) != 0.0)
                    Eint_New = Eint_Old - (Temp1 - Temp0) / ((Temp2 - Temp1) / (0.001 * Eint_Old));
                else
                    Eint_New = Eint_Old;

                double Epsilon = FABS(Eint_New - Eint_Old) / Eint_Old;
                Eint_Old = Eint_New;

                if (FABS(Epsilon) < 1e-4)
                    break;
                else if (ii == 999)
                {
					fprintf(stdout, "Newton for Eint(Dens,Temp) did not converge\n");
                }
            }
            DensTemp2Eint_Table[it * nRho + ir] = Eint_Old;
        }
    }
    fprintf(stdout, "%s ... done\n", __FUNCTION__);
}

void Construct_DensPres2Eint_Table()
{
    fprintf(stdout, "%s ... \n", __FUNCTION__);
    double Pres1, Pres1_CGS, Pres2, Pres2_CGS, Eint_New;

    DensPres2Eint_Table = (double *)malloc(sizeof(double) * nRho * nPres);
    memset(DensPres2Eint_Table, 0, sizeof(double) * nRho * nPres);

    for (int ir = 1; ir < nRho - 1; ir++)
    {
        for (int ip = 0; ip < nPres; ip++)
        {
            double Rho0 = POW(10.0, Dens_Table[1 * nRho + ir]);
            double Pres0 = POW(10.0, ( LOG10(PresMin) + ip * dPres)) * Rho0;
            double Eint_Old = Pres0 / (GAMMA - 1.0);
            
            //if (ip > 0) Eint_Old = MAX(Eint_Old, DensPres2Eint_Table[(ip - 1) * nRho + ir]);
            if (ip > 0) 
            {
                if (DensPres2Eint_Table[(ip - 1) * nRho + ir] != 0.0 )
                {
                    Eint_Old = DensPres2Eint_Table[(ip - 1) * nRho + ir];
                }
            }
            
            for (int ii = 0; ii < 1000; ii++)
            {   
                Pres1 = EoS_DensEint2Pres(Rho0/UNIT_D, Eint_Old/(UNIT_D * SQR(UNIT_V)));
                if (Pres1 == 0.0)
                {
                    Eint_Old = 0.0;
                    break;
                }

                Pres2 = EoS_DensEint2Pres(Rho0/UNIT_D, Eint_Old * 1.001/(UNIT_D * SQR(UNIT_V)));
                if (Pres2 == 0.0)
                {
                    Eint_Old = 0.0;
                    break;
                }
                
                if (FABS(Pres2 - Pres1) != 0.0)
                {
                    Pres1_CGS = Pres1 * (UNIT_D * SQR(UNIT_V));
                    Pres2_CGS = Pres2 * (UNIT_D * SQR(UNIT_V));
                    Eint_New = Eint_Old - (Pres1_CGS - Pres0 ) / ( (Pres2_CGS - Pres1_CGS ) / (0.001 * Eint_Old ));
                }
                else
                {
                    Eint_New = Eint_Old;
                }
                
                double Epsilon = FABS(Eint_New - Eint_Old) / Eint_Old;
                Eint_Old = Eint_New;

                if (FABS(Epsilon) < 1e-4)
                    break;
                else if (ii == 999)
                {
                    fprintf(stdout, "Newton for Eint(Dens,Pres) did not converge\n");
                }
            }

            DensPres2Eint_Table[ip * nRho + ir] = Eint_Old;
        }
    }

	fprintf(stdout, "%s ... done\n", __FUNCTION__);
}

