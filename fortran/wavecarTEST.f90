! -------1---------2---------3---------4---------5---------6---------7 !
! -------------------------------------------------------------------- !
! Hello!
! This bastardized script originated from 'wavetransSpinor.f90' code.
! It basically serves to write out one wavefunction in real-space,
! that is specified to by you!
!
! For example, (after you've compiled 'parseWAVECAR.exe'):
!    ./parseWAVECAR -f WAVECAR -k 80 -b 45
!
! (Defaults are -f WAVECAR -k 1 -b 1)
! This would read the WAVECAR file within the directory, then output
! the file WAVECAR_k80_b45.dat, with data preceding a 25-line header.
!
!
! and as always...
!                    ...stay sexy.
!
! Love,
!   Tyler
!
! -------------------------------------------------------------------- !
!  Jonathan T Reichanadter
!  jtreichanadter at berkeley.edu
!  Neaton Group - LBNL and UC Berkeley
!  September 15th, 2020
! -------------------------------------------------------------------- !
PROGRAM wavecarTEST



! --- 0 --- load modules
USE vtylib__WAVECAR
USE vtylib__POSCAR

! --- 1 --- workspace
IMPLICIT none
TYPE(userArgumentStruct) :: user
TYPE(wavecarFile) :: wavecar
TYPE(waveFunction) :: psi , psiTMP
TYPE(poscarFile) :: poscar
wavecar%unitNO = 69
poscar%unitNO = 12

! --- 2 --- parse data, and check wavecar
CALL parse( user )
CALL loadWAVECAR( user , wavecar )

! test poscar
poscar%fileName = "JOB__mos2mono_ncl/POSCAR"
CALL readPOSCAR( poscar )
CALL printPOSCAR( poscar )
CALL printWAVECAR( wavecar )

! --- 3 --- filter for exclusively noncollinear wavecars 
if (wavecar%magState.ne.'noncollinear') then
	write(0,*) "ERROR! This code requires a noncollinear WAVECAR"
	stop
endif


! --- 4 --- read seed wavefunction
CALL loadWAVEFUNCTION( wavecar , user%k , user%b , psi )



!do q = 0,(user%kMAX-user%kMIN)
!	do m = 0,(user%bMAX-user%bMIN)
!		CALL loadWAVEFUNCTION( wavecar , kMIN+q , bMIN+m , psi1 )
!			do j = 1,wavecarFile%npCAP
!
!			enddo
!	enddo
!enddo



















end PROGRAM
