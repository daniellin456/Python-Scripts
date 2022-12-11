#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "Create_Table.h"


int nRho;
int nEnergy;
int nTemp;
int nPres;

const double UNIT_L  = 3.08e18;
const double UNIT_T  = 1.9773e+15;
const double UNIT_V  = 3.08e18/1.9773e+15;
const double UNIT_D  = 3.8346e-24;
const double Const_kB = 1.38064852e-16;
const double Const_amu = 1.660539040e-24;
const double GAMMA = 1.666666667;
const double MOLECULAR_WEIGHT = 2.31;

double RhoMin;
double RhoMax;
double Emin;
double Emax;
double TempMin;
double TempMax;
double PresMin;
double PresMax;

double dRho;
double dEner;
double dTemp;
double dPres;
double yHe;

double *Dens_Table;
double *Ener_Table;
double *DensEint2Pres_Table;
double *DensEint2Temp_Table;
double *DensEint2Cs_Table;
double *DensPres2Eint_Table;
double *DensTemp2Eint_Table;

void Dump_File( char *Filename, double* table, const int nRho, const int nEnergy )
{
	FILE *fp = NULL;

	if ( (fp = fopen(Filename, "wb")) == NULL )
	{
		printf("Could not open <%s>.\n", Filename);
		exit(1);
	}

	fwrite(table, sizeof(double) * nRho * nEnergy, 1, fp);
	fclose(fp);
}

void Fill_Hole_in_Table(double *table, const int nRho, const int nEnergy)
{
	int ii, jj, kk, hh, ee, gg;
	for (int k = 0; k < 5; k++)
	{
		ii = 0;
		jj = 0;
		kk = 0;
		hh = 0;
		ee = 0;
		gg = 0;
		double ww, xx, yy, zz;

		for (int ir = 1; ir < nRho - 1; ir++)
		{
			for (int ie = 1; ie < nEnergy - 1; ie++)
			{
				if (table[ie * nRho + ir] == 0.0)
				{
					ii++;
					xx = table[(ie + 1) * nRho + ir] *
						 table[(ie - 1) * nRho + ir] *
						 table[ie * nRho + ir - 1] *
						 table[ie * nRho + ir + 1];

					yy = table[(ie + 1) * nRho + ir + 1] *
						 table[(ie - 1) * nRho + ir + 1] *
						 table[(ie - 1) * nRho + ir - 1] *
						 table[(ie + 1) * nRho + ir - 1];

					if (ie > 1 && ie < nEnergy - 2 && ir > 1 && ir < nRho - 2)
					{
						ww = table[(ie + 2) * nRho + ir] *
							 table[(ie - 2) * nRho + ir] *
							 table[ie * nRho + ir - 2] *
							 table[ie * nRho + ir + 2];
					}
					else
						ww = 0.0;

					if (ie > 2 && ie < nEnergy - 3 && ir > 2 && ir < nRho - 3)
					{
						zz = table[(ie + 3) * nRho + ir + 3] *
							 table[(ie - 3) * nRho + ir - 3] *
							 table[(ie + 3) * nRho + ir - 3] *
							 table[(ie - 3) * nRho + ir + 3];
					}
					else
						zz = 0.0;

					if (xx != 0.0)
					{
						jj++;
						table[ie * nRho + ir] = 0.25 * (table[(ie + 1) * nRho + ir] +
														table[(ie - 1) * nRho + ir] +
														table[ie * nRho + ir - 1] +
														table[ie * nRho + ir + 1]);
					}
					else if (yy != 0.0 && k >= 0)
					{
						kk++;
						table[ie * nRho + ir] = 0.25 * (table[(ie + 1) * nRho + ir + 1] +
														table[(ie - 1) * nRho + ir + 1] +
														table[(ie + 1) * nRho + ir - 1] +
														table[(ie - 1) * nRho + ir - 1]);
					}
					else if (ww != 0 && k > 0)
					{
						ee++;
						table[ie * nRho + ir] = 0.25 * (table[(ie + 2) * nRho + ir] +
														table[(ie - 2) * nRho + ir] +
														table[ie * nRho + ir - 2] +
														table[ie * nRho + ir + 2]);
					}
					else if (zz != 0 && k > 1)
					{
						hh++;
						table[ie * nRho + ir] = 0.25 * (table[(ie + 3) * nRho + ir + 3] +
														table[(ie - 3) * nRho + ir + 3] +
														table[(ie + 3) * nRho + ir - 3] +
														table[(ie - 3) * nRho + ir - 3]);
					}
					else
						gg++;
				}
			}
		}
		printf("Fill the hole of Table, ii = %5d, jj = %5d, kk = %5d, ee = %5d, hh = %5d, gg = %5d, iter = %5d\n",
			   ii, jj, kk, ee, hh, gg, k);
	}
}

int main(int argc, char *argv[])
{
	FILE *fp = NULL;
	char *EoS_Table_Filename =						(char *)malloc(50 * sizeof(char));
	char *Dens_Table_Filename =						(char *)malloc(50 * sizeof(char));
	char *Ener_Table_Filename =						(char *)malloc(50 * sizeof(char));
	char *DensEint2Pres_Table_Filename = 			(char *)malloc(50 * sizeof(char));
	char *DensEint2Temp_Table_Filename = 			(char *)malloc(50 * sizeof(char));
	char *DensEint2Cs_Table_Filename =   			(char *)malloc(50 * sizeof(char));
	char *Fill_DensEint2Pres_Table_Filename = 		(char *)malloc(50 * sizeof(char));
	char *Fill_DensEint2Temp_Table_Filename = 		(char *)malloc(50 * sizeof(char));
	char *Fill_DensEint2Cs_Table_Filename =   		(char *)malloc(50 * sizeof(char));
	char *DensPres2Eint_Table_Filename   =          (char *)malloc(50 * sizeof(char));
	char *DensTemp2Eint_Table_Filename   =          (char *)malloc(50 * sizeof(char));

	strcpy( EoS_Table_Filename,						"../resources/tab_eos.dat"					);
	strcpy( Dens_Table_Filename,					"../resources/Dens_Table.bin"				);
	strcpy( Ener_Table_Filename,					"../resources/Ener_Table.bin"				);
	strcpy( DensEint2Pres_Table_Filename, 		    "../resources/DensEint2Pres_Table.bin"		);
	strcpy( DensEint2Temp_Table_Filename, 		    "../resources/DensEint2Temp_Table.bin"		);
	strcpy( DensEint2Cs_Table_Filename,   		    "../resources/DensEint2Cs_Table.bin"		);
	strcpy( Fill_DensEint2Pres_Table_Filename, 	    "../resources/Fill_DensEint2Pres_Table.bin"	);
	strcpy( Fill_DensEint2Temp_Table_Filename, 	    "../resources/Fill_DensEint2Temp_Table.bin"	);
	strcpy( Fill_DensEint2Cs_Table_Filename, 	    "../resources/Fill_DensEint2Cs_Table.bin"	);
	strcpy( DensPres2Eint_Table_Filename,			"../resources/DensPres2Eint_Table.bin"		);
	strcpy( DensTemp2Eint_Table_Filename,			"../resources/DensTemp2Eint_Table.bin"		);
	
	if ( ( fp = fopen(EoS_Table_Filename, "rb") ) == NULL)
	{
		printf("Could not open <%s>.\n", EoS_Table_Filename);
		return 1;
	}

	fseek(fp, 0, SEEK_SET); // Move cursor to begin of file
	fseek(fp, 4, SEEK_CUR); // Skip the blank

	fread(&nRho, sizeof(int), 1, fp);
	fread(&nEnergy, sizeof(int), 1, fp);

	fseek(fp, 8, SEEK_CUR); // Skip the blank
	fread(&RhoMin, sizeof(double), 1, fp);
	fread(&RhoMax, sizeof(double), 1, fp);
	fread(&Emin, sizeof(double), 1, fp);
	fread(&Emax, sizeof(double), 1, fp);
	fread(&yHe, sizeof(double), 1, fp);

	fseek(fp, 8, SEEK_CUR); // Skip the blank
	Dens_Table = (double *)malloc(sizeof(double) * nRho * nEnergy);
	fread(Dens_Table, sizeof(double) * nRho * nEnergy, 1, fp);

	fseek(fp, 8, SEEK_CUR); // Skip the blank
	Ener_Table = (double *)malloc(sizeof(double) * nRho * nEnergy);
	fread(Ener_Table, sizeof(double) * nRho * nEnergy, 1, fp);

	fseek(fp, 8, SEEK_CUR); // Skip the blank
	DensEint2Temp_Table = (double *)malloc(sizeof(double) * nRho * nEnergy);
	fread(DensEint2Temp_Table, sizeof(double) * nRho * nEnergy, 1, fp);

	fseek(fp, 8, SEEK_CUR); // Skip the blank
	DensEint2Pres_Table = (double *)malloc(sizeof(double) * nRho * nEnergy);
	fread(DensEint2Pres_Table, sizeof(double) * nRho * nEnergy, 1, fp);

	fseek(fp, 8, SEEK_CUR); // Skip the blank
	fseek(fp, sizeof(double) * nRho * nEnergy, SEEK_CUR ); // Skip S

	fseek(fp, 8, SEEK_CUR); // Skip the blank
	DensEint2Cs_Table = (double *)malloc(sizeof(double) * nRho * nEnergy);
	fread(DensEint2Cs_Table, sizeof(double) * nRho * nEnergy, 1, fp);

	fclose(fp); // fopen("tab_eos.dat")

	for( int i = 0 ; i < nRho * nEnergy; i++ )
	{
		Dens_Table[i] = log10(Dens_Table[i]);
		Ener_Table[i] = log10(Ener_Table[i]);
	}

	printf("Ener_Table[0] = %20.14e\n", Ener_Table[0]);
	printf("Ener_Table[623] = %20.14e\n", Ener_Table[623]);
	printf("Ener_Table[624] = %20.14e\n", Ener_Table[624]);
	printf("Ener_Table[625] = %20.14e\n", Ener_Table[625*999-1]);
	printf("Ener_Table[1249] = %20.14e\n", Ener_Table[625*1000-1]);

	nTemp = nEnergy;
	nPres = nEnergy;
    TempMax = 1e5;
    TempMin = 3.0;
    PresMax = 1e14;
    PresMin = 1e8;

	dRho  = ( RhoMax - RhoMin ) / nRho;
	dEner = ( Emax - Emin ) / nEnergy ;
	dTemp = ( log10(TempMax) - log10(TempMin) ) / nTemp;
	dPres = ( log10(PresMax) - log10(PresMin) ) / nPres;

	Dump_File( Dens_Table_Filename, Dens_Table, nRho, nEnergy );
	Dump_File( Ener_Table_Filename, Ener_Table, nRho, nEnergy );
	Dump_File( DensEint2Temp_Table_Filename, DensEint2Temp_Table, nRho, nEnergy );
	Dump_File( DensEint2Pres_Table_Filename, DensEint2Pres_Table, nRho, nEnergy );
	Dump_File( DensEint2Cs_Table_Filename,   DensEint2Cs_Table,   nRho, nEnergy );

	Fill_Hole_in_Table( DensEint2Pres_Table, nRho, nEnergy );
	Fill_Hole_in_Table( DensEint2Temp_Table, nRho, nEnergy );
	Fill_Hole_in_Table( DensEint2Cs_Table,   nRho, nEnergy );

	Dump_File( Fill_DensEint2Pres_Table_Filename, DensEint2Pres_Table, nRho, nEnergy );
	Dump_File( Fill_DensEint2Temp_Table_Filename, DensEint2Temp_Table, nRho, nEnergy );
	Dump_File( Fill_DensEint2Cs_Table_Filename,   DensEint2Cs_Table,   nRho, nEnergy );

	Construct_DensPres2Eint_Table();
	Construct_DensTemp2Eint_Table();

	Dump_File( DensPres2Eint_Table_Filename, DensPres2Eint_Table, nRho, nEnergy );
	Dump_File( DensTemp2Eint_Table_Filename, DensTemp2Eint_Table, nRho, nEnergy );
	
	return 0;
}
