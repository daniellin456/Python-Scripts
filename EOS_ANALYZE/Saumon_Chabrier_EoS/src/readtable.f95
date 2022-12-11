PROGRAM readTable

USE CONSTANTS

IMPLICIT NONE

REAL( KIND = 8 ) :: dTemp1, Rho0, T0, Eint_Old, Eint_New, Temp_New, Temp_New2, epsilon
REAL( KIND = 8 ) :: ww, xx, yy, zz
INTEGER :: ii, jj, kk, hh, ee, gg, i, j, k, ir, it, ht

OPEN( 10, FILE = "../resources/tab_eos.dat", STATUS = "UNKNOWN", FORM = "UNFORMATTED" )
READ( 10 ) nRho, nEnergy
READ( 10 ) RhoMin, RhoMax, Emin, Emax, yHe

ALLOCATE( Dens_Table( nRho, nEnergy ) )
ALLOCATE( Ener_Table( nRho, nEnergy ) )
ALLOCATE( DensEint2Pres_Table( nRho, nEnergy ) )
ALLOCATE( DensEint2Temp_Table( nRho, nEnergy ) )
ALLOCATE( S_Table( nRho, nEnergy ) )
ALLOCATE( DensEint2Cs_Table( nRho, nEnergy ) )

READ( 10 ) Dens_Table
READ( 10 ) Ener_Table
READ( 10 ) DensEint2Temp_Table
READ( 10 ) DensEint2Pres_Table
READ( 10 ) S_Table
READ( 10 ) DensEint2Cs_Table

CLOSE( 10 )

! Fill the hole DensEint2Pres_Table
DO k = 1, 5
    ii = 0
    jj = 0 
    kk = 0
    hh = 0
    ee = 0
    gg = 0
    
    DO i = 2, nRho - 1
        DO j = 2, nEnergy - 1
            IF( DensEint2Pres_Table( i, j ) .EQ. 0.0d0 ) THEN
                ii = ii + 1
                xx = DensEint2Pres_Table( i, j + 1 ) * &
                     DensEint2Pres_Table( i, j - 1 ) * &
                     DensEint2Pres_Table( i - 1, j ) * &
                     DensEint2Pres_Table( i + 1, j )
                yy = DensEint2Pres_Table( i + 1, j + 1 ) * &
                     DensEint2Pres_Table( i + 1, j - 1 ) * &
                     DensEint2Pres_Table( i - 1, j - 1 ) * &
                     DensEint2Pres_Table( i - 1, j + 1 )

                IF( j > 2 .AND. j < nEnergy - 1 .AND. & 
                    i > 2 .AND. i < nRho - 1 ) THEN
                    ww = DensEint2Pres_Table( i, j + 2 ) * &
                         DensEint2Pres_Table( i, j - 2 ) * &
                         DensEint2Pres_Table( i - 2, j ) * &
                         DensEint2Pres_Table( i + 2, j )
                ELSE
                    ww = 0.0
                END IF

                IF( j > 3 .AND. j < nEnergy -2 .AND. &
                    i > 3 .AND. i < nRho - 2 ) THEN
                    zz = DensEint2Pres_Table( i + 3, j + 3 ) * &
                         DensEint2Pres_Table( i - 3, j - 3 ) * &
                         DensEint2Pres_Table( i - 3, j + 3 ) * &
                         DensEint2Pres_Table( i + 3, j - 3 ) 
                ELSE
                    zz = 0.0
                END IF

                IF( xx .NE. 0 ) THEN
                    jj = jj + 1
                    DensEint2Pres_Table( i, j ) = 0.25 * ( DensEint2Pres_Table( i, j + 1 ) + &
                                                  DensEint2Pres_Table( i, j - 1 ) + &
                                                  DensEint2Pres_Table( i - 1, j ) + &
                                                  DensEint2Pres_Table( i + 1, j ) )
                ELSE IF( yy .NE. 0 .AND. k > 0 ) THEN
                    kk = kk + 1
                    DensEint2Pres_Table( i, j ) = 0.25 * ( DensEint2Pres_Table( i + 1, j + 1 ) + &
                                                  DensEint2Pres_Table( i + 1, j - 1 ) + &
                                                  DensEint2Pres_Table( i - 1, j + 1 ) + &
                                                  DensEint2Pres_Table( i - 1, j - 1 ) )

                ELSE IF( ww .NE. 0 .AND. k > 1 ) THEN
                    ee = ee + 1
                    DensEint2Pres_Table( i, j ) = 0.25 * ( DensEint2Pres_Table( i, j + 2 ) + &
                                                  DensEint2Pres_Table( i, j - 2 ) + &
                                                  DensEint2Pres_Table( i - 2, j ) + &
                                                  DensEint2Pres_Table( i + 2, j ) )

                ELSE IF( zz .NE. 0 .AND. k > 2 ) THEN
                    hh = hh + 1
                    DensEint2Pres_Table( i, j ) = 0.25 * ( DensEint2Pres_Table( i + 3, j + 3 ) + &
                                                  DensEint2Pres_Table( i + 3, j - 3 ) + &
                                                  DensEint2Pres_Table( i - 3, j + 3 ) + &
                                                  DensEint2Pres_Table( i - 3, j - 3 ) )
                ELSE
                    gg = gg + 1
                END IF
            END IF
        END DO
    END DO
    WRITE( *, * ) "on bouche les trous DensEint2Pres_Table", ii,jj,kk,ee,hh,gg, "iter", k
END DO

! Fill the hole DensEint2Cs_Table
DO k = 1, 5
    ii = 0
    jj = 0 
    kk = 0
    hh = 0
    ee = 0
    gg = 0
    
    DO i = 2, nRho - 1
        DO j = 2, nEnergy - 1
            IF( DensEint2Cs_Table( i, j ) .EQ. 0.0d0 ) THEN
                ii = ii + 1
                xx = DensEint2Cs_Table( i, j + 1 ) * &
                     DensEint2Cs_Table( i, j - 1 ) * &
                     DensEint2Cs_Table( i - 1, j ) * &
                     DensEint2Cs_Table( i + 1, j )
                yy = DensEint2Cs_Table( i + 1, j + 1 ) * &
                     DensEint2Cs_Table( i + 1, j - 1 ) * &
                     DensEint2Cs_Table( i - 1, j - 1 ) * &
                     DensEint2Cs_Table( i - 1, j + 1 )

                IF( j > 2 .AND. j < nEnergy - 1 .AND. & 
                    i > 2 .AND. i < nRho - 1 ) THEN
                    ww = DensEint2Cs_Table( i, j + 2 ) * &
                         DensEint2Cs_Table( i, j - 2 ) * &
                         DensEint2Cs_Table( i - 2, j ) * &
                         DensEint2Cs_Table( i + 2, j )
                ELSE
                    ww = 0.0
                END IF

                IF( j > 3 .AND. j < nEnergy -2 .AND. &
                    i > 3 .AND. i < nRho - 2 ) THEN
                    zz = DensEint2Cs_Table( i + 3, j + 3 ) * &
                         DensEint2Cs_Table( i - 3, j - 3 ) * &
                         DensEint2Cs_Table( i - 3, j + 3 ) * &
                         DensEint2Cs_Table( i + 3, j - 3 ) 
                ELSE
                    zz = 0.0
                END IF

                IF( xx .NE. 0 ) THEN
                    jj = jj + 1
                    DensEint2Cs_Table( i, j ) = 0.25 * ( DensEint2Cs_Table( i, j + 1 ) + &
                                                  DensEint2Cs_Table( i, j - 1 ) + &
                                                  DensEint2Cs_Table( i - 1, j ) + &
                                                  DensEint2Cs_Table( i + 1, j ) )
                ELSE IF( yy .NE. 0 .AND. k > 0 ) THEN
                    kk = kk + 1
                    DensEint2Cs_Table( i, j ) = 0.25 * ( DensEint2Cs_Table( i + 1, j + 1 ) + &
                                                  DensEint2Cs_Table( i + 1, j - 1 ) + &
                                                  DensEint2Cs_Table( i - 1, j + 1 ) + &
                                                  DensEint2Cs_Table( i - 1, j - 1 ) )

                ELSE IF( ww .NE. 0 .AND. k > 1 ) THEN
                    ee = ee + 1
                    DensEint2Cs_Table( i, j ) = 0.25 * ( DensEint2Cs_Table( i, j + 2 ) + &
                                                  DensEint2Cs_Table( i, j - 2 ) + &
                                                  DensEint2Cs_Table( i - 2, j ) + &
                                                  DensEint2Cs_Table( i + 2, j ) )

                ELSE IF( zz .NE. 0 .AND. k > 2 ) THEN
                    hh = hh + 1
                    DensEint2Cs_Table( i, j ) = 0.25 * ( DensEint2Cs_Table( i + 3, j + 3 ) + &
                                                  DensEint2Cs_Table( i + 3, j - 3 ) + &
                                                  DensEint2Cs_Table( i - 3, j + 3 ) + &
                                                  DensEint2Cs_Table( i - 3, j - 3 ) )
                ELSE
                    gg = gg + 1
                END IF
            END IF
        END DO
    END DO
    WRITE( *, * ) "on bouche les trous DensEint2Cs_Table", ii,jj,kk,ee,hh,gg, "iter", k
END DO

! Fill the hole DensEint2Temp_Table
DO k = 1, 5
    ii = 0
    jj = 0 
    kk = 0
    hh = 0
    ee = 0
    gg = 0
    
    DO i = 2, nRho - 1
        DO j = 2, nEnergy - 1
            IF( DensEint2Temp_Table( i, j ) .EQ. 0.0d0 ) THEN
                ii = ii + 1
                xx = DensEint2Temp_Table( i, j + 1 ) * &
                     DensEint2Temp_Table( i, j - 1 ) * &
                     DensEint2Temp_Table( i - 1, j ) * &
                     DensEint2Temp_Table( i + 1, j )
                yy = DensEint2Temp_Table( i + 1, j + 1 ) * &
                     DensEint2Temp_Table( i + 1, j - 1 ) * &
                     DensEint2Temp_Table( i - 1, j - 1 ) * &
                     DensEint2Temp_Table( i - 1, j + 1 )

                IF( j > 2 .AND. j < nEnergy - 1 .AND. & 
                    i > 2 .AND. i < nRho - 1 ) THEN
                    ww = DensEint2Temp_Table( i, j + 2 ) * &
                         DensEint2Temp_Table( i, j - 2 ) * &
                         DensEint2Temp_Table( i - 2, j ) * &
                         DensEint2Temp_Table( i + 2, j )
                ELSE
                    ww = 0.0
                END IF

                IF( j > 3 .AND. j < nEnergy -2 .AND. &
                    i > 3 .AND. i < nRho - 2 ) THEN
                    zz = DensEint2Temp_Table( i + 3, j + 3 ) * &
                         DensEint2Temp_Table( i - 3, j - 3 ) * &
                         DensEint2Temp_Table( i - 3, j + 3 ) * &
                         DensEint2Temp_Table( i + 3, j - 3 ) 
                ELSE
                    zz = 0.0
                END IF

                IF( xx .NE. 0 ) THEN
                    jj = jj + 1
                    DensEint2Temp_Table( i, j ) = 0.25 * ( DensEint2Temp_Table( i, j + 1 ) + &
                                                  DensEint2Temp_Table( i, j - 1 ) + &
                                                  DensEint2Temp_Table( i - 1, j ) + &
                                                  DensEint2Temp_Table( i + 1, j ) )
                ELSE IF( yy .NE. 0 .AND. k > 0 ) THEN
                    kk = kk + 1
                    DensEint2Temp_Table( i, j ) = 0.25 * ( DensEint2Temp_Table( i + 1, j + 1 ) + &
                                                  DensEint2Temp_Table( i + 1, j - 1 ) + &
                                                  DensEint2Temp_Table( i - 1, j + 1 ) + &
                                                  DensEint2Temp_Table( i - 1, j - 1 ) )

                ELSE IF( ww .NE. 0 .AND. k > 1 ) THEN
                    ee = ee + 1
                    DensEint2Temp_Table( i, j ) = 0.25 * ( DensEint2Temp_Table( i, j + 2 ) + &
                                                  DensEint2Temp_Table( i, j - 2 ) + &
                                                  DensEint2Temp_Table( i - 2, j ) + &
                                                  DensEint2Temp_Table( i + 2, j ) )

                ELSE IF( zz .NE. 0 .AND. k > 2 ) THEN
                    hh = hh + 1
                    DensEint2Temp_Table( i, j ) = 0.25 * ( DensEint2Temp_Table( i + 3, j + 3 ) + &
                                                  DensEint2Temp_Table( i + 3, j - 3 ) + &
                                                  DensEint2Temp_Table( i - 3, j + 3 ) + &
                                                  DensEint2Temp_Table( i - 3, j - 3 ) )
                ELSE
                    gg = gg + 1
                END IF
            END IF
        END DO
    END DO
    WRITE( *, * ) "on bouche les trous DensEint2Temp_Table", ii,jj,kk,ee,hh,gg, "iter", k
END DO

nTemp = nEnergy
ALLOCATE( DensTemp2Eint_Table( nRho, nTemp ) )
dTemp1 = ( log10(Tmax) - log10(Tmin) ) / nTemp
DensTemp2Eint_Table(:,:)=0.0d0

DO ir = 2, nRho - 1
    DO it = 1, nTemp
        
        Rho0 = ( 10.0 ** Dens_Table(ir,1) )
        T0 = 10.0 ** (log10(Tmin) + ( it - 1.0d0 ) * dtemp1 )
        Eint_Old = Rho0 * kb * T0 / ( mu_gas * mh * ( gamma - 1.0d0 ) )

        IF ( it > 1 ) THEN
            Eint_Old = max( Rho0 * kB * T0 / ( mu_gas * mH * ( gamma - 1.0d0 ) ), Ener_Table(ir,it-1) )
        END IF
      
        epsilon = 1.0d0
      
        DO ii = 1, 1000
            
            CALL temperature_eos( Rho0 / scale_d, Eint_Old / ( scale_d * scale_v**2 ), Temp_New, ht )
            IF ( ht == 1 ) THEN
                Eint_Old = 0.d0
                EXIT
            END IF

            CALL temperature_eos( Rho0 / scale_d, Eint_Old * 1.001 / ( scale_d * scale_v**2 ), Temp_New2, ht )
            IF ( ht == 1 ) THEN
                Eint_Old = 0.d0
                EXIT
            END IF

            IF( abs( Temp_New2 - Temp_New ) .NE. 0 ) THEN
                Eint_New = Eint_Old - ( Temp_New - T0 ) / ( ( Temp_New2 - Temp_New ) / ( 0.001 * Eint_Old ) )
            ELSE
                Eint_New = Eint_Old
            END IF

            epsilon = abs(Eint_New - Eint_Old)/Eint_Old
            Eint_Old = Eint_New
            
            IF ( abs(epsilon) .LT. 1.d-4) THEN
                EXIT
            ELSE IF  ( ii == 1000 ) THEN
                PRINT *, "newton for e(rho,T) did not converge at ", log10(Rho0), log10(T0)
            END IF

        END DO 
        DensTemp2Eint_Table(ir,it) = Eint_Old 
    END DO
END DO

OPEN( 10, FILE = "../resources/DensTemp2Eint_Table_Fortran.bin", STATUS = "UNKNOWN", FORM = "UNFORMATTED" )
WRITE( 10 ) DensTemp2Eint_Table
CLOSE( 10 )

END PROGRAM readTable


SUBROUTINE temperature_eos( Dens, Eint, Temp, ht )
    
    USE CONSTANTS

    IMPLICIT NONE

    REAL( KIND = 8 ), intent(in) :: Eint, Dens
    REAL( KIND = 8 ), intent(out):: Temp
    REAL( KIND = 8 ) :: Eint_CGS, Dens_CGS
    REAL( KIND = 8 ) :: le, lr
    REAL( KIND = 8 ) :: dd1, dd2, de1, de2
    REAL( KIND = 8 ) :: xx, dDens, dEner
    INTEGER :: ht
    INTEGER :: ir,ie

    IF (Eint == 0.d0) THEN
        Temp = 0.d0
    ELSE
        ht=0
        Dens_CGS = Dens * scale_d             
        Eint_CGS = Eint * scale_d * scale_v**2
   
        dDens = ( RhoMax - RhoMin ) / float( nRho )
        lr = 0.5d0 + ( log10(Dens_CGS) - RhoMin ) / dDens
   
        IF ( lr .GE. nRho ) THEN
           WRITE(*,*)'pb 1'
           WRITE(*,*) log10( Dens_CGS )
           STOP
        END IF

        ir = floor(lr)
   
        dEner = ( Emax  -   Emin ) / float( nEnergy )
        le = 0.5d0 + (log10(Eint_CGS) - Emin - log10(Dens_CGS) )/ dEner
   
        IF ((le .GE. nEnergy) ) THEN
           WRITE(*,*)'pb 2'
           STOP
        END IF
   
        ie = floor(le)
        IF  ( ir < 1 .OR. ie < 1 ) THEN
            WRITE(*,*) 'inter_tp hors limite ir,ie,rho,enint = ', ir, ie, Dens, Eint 
            ir=1.0d0
            ie=1.0d0
            STOP
        END IF
   
        dd1 = lr - float(ir)
        dd2 = le - float(ie)
   
        de1 = 1.0d0 - dd1
        de2 = 1.0d0 - dd2
   
        Temp = 0.d0
   
        Temp = Temp + de1 * de2 * DensEint2Temp_Table(ir  ,ie  )
        Temp = Temp + dd1 * de2 * DensEint2Temp_Table(ir+1,ie  )
        Temp = Temp + de1 * dd2 * DensEint2Temp_Table(ir  ,ie+1)
        Temp = Temp + dd1 * dd2 * DensEint2Temp_Table(ir+1,ie+1)
   
        Temp = Temp ! give T in K
   
        xx =  DensEint2Temp_Table(ir,ie)*DensEint2Temp_Table(ir+1,ie)*DensEint2Temp_Table(ir,ie+1)*DensEint2Temp_Table(ir+1,ie+1) 
   
        IF (xx .EQ. 0.0d0 ) THEN
           ht=1
        END IF
    END IF
END SUBROUTINE