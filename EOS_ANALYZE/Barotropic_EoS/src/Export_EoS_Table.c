#include <stdio.h>
#include <stdlib.h>

int main( int argc, char *argv[] )
{
	FILE *fp = NULL;
	char filename[50] = "barotropic_eos.dat";
	char DensEoS[50] = "DensEoS.csv";
	char TempEoS[50] = "TempEoS.csv";
	int nRho;
	double RhoMin, RhoMax, dRho;
	double *Dens_Table, *Temp_Table;

	fp = fopen( filename, "rb" );
	if ( fp == NULL )
	{
		printf("Could not open <%s>.\n", filename );
		return 1;
	}
	fscanf( fp, "\t%d\t%lf\t%lf\t%lf", &nRho, &RhoMin, &RhoMax, &dRho );
	fprintf( stdout, "%d %lf %lf %lf\n", nRho, RhoMin, RhoMax, dRho);

	Dens_Table = ( double * ) malloc ( sizeof( double ) * nRho );
	Temp_Table = ( double * ) malloc ( sizeof( double ) * nRho );
	for ( int i = 0; i < nRho; i++ )
	{
		fscanf(fp, "\t%lf\t%lf", &Dens_Table[i], &Temp_Table[i] );
	}

	/*
	for ( int i = 0; i < nRho; i++ )
	{
		fprintf( stdout, "%20.14e\t%20.14e\n", Dens_Table[i], Temp_Table[i] );
	}
	*/

	fclose( fp );
	
	fp = fopen( DensEoS, "w" );
	if ( fp == NULL )
	{
		printf("Could not open <%s>.\n", DensEoS );
		return 1;
	}
	for ( int i = 0; i < nRho; i++ ) fprintf( fp, "%20.14e\n", Dens_Table[i] );
	fclose( fp );
	
	fp = fopen( TempEoS, "w" );
	if ( fp == NULL )
	{
		printf("Could not open <%s>.\n", TempEoS );
		return 1;
	}
	for ( int i = 0; i < nRho; i++ ) fprintf( fp, "%20.14e\n", Temp_Table[i] );
	fclose( fp );

	return 0;
}
