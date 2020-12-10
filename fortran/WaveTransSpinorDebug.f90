!!$************************* WaveTransSpinor *******************************
!!$
!!$   input the WAVECAR file in binary format from VASP, and write
!!$   selected real space wavefunction in a3 direction to standard output
!!$   Plane wave coefficients are written to GCOEFF.txt
!!$
!!$   Compile with gfortran or ifort. Flag "-assume byterecl" is required
!!$   for ifort.
!!$
!!$   version 2.0 - Jul 11, 2012 - R. M. Feenstra and M. Widom, 
!!$                                using updated 'c' value
!!$   version 2.1 - Sep 17, 2014 - changed estimator for max. no. of
!!$                                plane waves
!!$   version 2.2 - Nov 12, 2014 - output (n1,n2,n3) indices in both first half
!!$                                and second half of GCOEFF.txt output file
!!$
!!$   options are -f filename -k k-point -b band -x coord -y coord
!!$   defaults are -f WAVECAR -k 1 -b 1 -x 0 -y 0
!!$   coordinates are direct coordinates



!!$* ----------------- initialize workspace variables ----------------------
!!$* allow implicit variable type-casting
implicit real*8 (a-h, o-z)
!!$* initialize dynamically allocated memory for data
complex*8, allocatable :: coeff(:)
complex*16, allocatable :: cener(:)
real*8, allocatable :: occ(:)
integer, allocatable :: igall(:,:)
!!$* initialize vector variables
dimension a1(3),a2(3),a3(3),b1(3),b2(3),b3(3),a2xa3(3),sumkg(3),vtmp(3)
dimension wk(3),xyz(3),wkpg(3),ig(3)
!!$* initialize intermediate variables
complex*16 csum1,csum2
integer kpoint,band
!!$* reserve variable for input file name
character*75 filename



!!$* -------------------------- set constants ------------------------------
!!$*   constant 'c' below is 2m/hbar**2 in units of 1/eV Ang^2 (value is
!!$*   adjusted in final decimal places to agree with VASP value; program
!!$*   checks for discrepancy of any results between this and VASP values)
data c/0.262465831d0/ 
pi=4.*atan(1.)
write(0,*) ''//achar(27)//'[0;38;5;021m' , '[0] --- Constants', ''//achar(27)//'[0m'
write(0,*) ''//achar(27)//'[0;38;5;033m' , 'c = ',c , ''//achar(27)//'[0m'
write(0,*) ' '




!!$ parse arguments
call parse(filename,kpoint,band,x,y)
xyz(1)=x
xyz(2)=y



!!$*  OPEN data and check record length
nrecl=24  ! "record length" for byte reading
open(unit=10,file=filename,access='direct',recl=nrecl,iostat=iost,status='old')
!!$* give up if we observe ANY sort of file opening error
if (iost.ne.0) write(6,*) 'open error - iostat =',iost            



!!$* READ out header variables, then CLOSE
read(unit=10,rec=1) xnrecl , xnspin , xnprec
close(unit=10)
nrecl=nint(xnrecl)
nspin=nint(xnspin)
nprec=nint(xnprec)

!!$ ---------- PRINT STATEMENTS -------------
write(0,*) ''//achar(27)//'[0;38;5;021m' , '[1] --- Precursors', ''//achar(27)//'[0m'
write(0,*) ''//achar(27)//'[0;38;5;021m' , '    read(unit=10,rec=1) xnrecl , xnspin , xnprec', ''//achar(27)//'[0m'
write(0,*) ''//achar(27)//'[0;38;5;033m' , 'nrecl = ',nrecl , ''//achar(27)//'[0m'
write(0,*) ''//achar(27)//'[0;38;5;033m' , 'nspin = ',nspin , ''//achar(27)//'[0m'
write(0,*) ''//achar(27)//'[0;38;5;033m' , 'nprec = ',nprec , ''//achar(27)//'[0m'
write(0,*) ' '










!! check data precision/format
if(nprec.eq.45210) then
   write(0,*) '*** error - WAVECAR_double requires complex*16'
   stop
endif
if(nspin.eq.2) then
   write(0,*) '*** error - Not a spinor WAVECAR. ISPIN =',nspin
   stop
endif










!$* REOPEN data file and CREATE 'gcoeff_mod.txt' file
!$* read out iostat ERRORS with our input
!!$*  10 - Read error on direct file.
!!$*  11 - Write error on direct file.
!!$*  37 - Dynamic memory allocation failure (out of memory).
!!$* See this link: www.ibm.com/support/knowledgecenter/en/SS3KZ4_9.0.0/com.ibm.xlf111.bg.doc/xlflr/iostatvalues.htm
open(unit=10,file=filename,access='direct',recl=nrecl, &
     iostat=iost,status='old')
if (iost.ne.0) write(6,*) 'open error - iostat =',iost
! open/create GCOEFF_MOD.txt file
open(unit=11,file='GCOEFF_MOD.txt')
read(unit=10,rec=2) xnwk,xnband,ecut,(a1(j),j=1,3),(a2(j),j=1,3),(a3(j),j=1,3)

!!$ ---------- PRINT STATEMENTS -------------
write(0,*) ''//achar(27)//'[0;38;5;021m' , '[1] --- Header Variables' , ''//achar(27)//'[0m'
write(0,*) ''//achar(27)//'[0;38;5;021m' , '    read(unit=10,rec=2) xnwk,xnband,ecut,' , &
'(a1(j),j=1,3),(a2(j),j=1,3),(a3(j),j=1,3)' , ''//achar(27)//'[0m'
write(0,*) ''//achar(27)//'[0;38;5;033m' , 'nwk = ' , nint(xnwk) , ''//achar(27)//'[0m'
write(0,*) ''//achar(27)//'[0;38;5;033m' , 'nband = ' , nint(xnband) , ''//achar(27)//'[0m'
write(0,*) ''//achar(27)//'[0;38;5;033m' , 'ecut = ' , nint(ecut) , ''//achar(27)//'[0m'
write(0,*) ''//achar(27)//'[0;38;5;033m' , 'a1(j) = ' , (a1(j),j=1,3) , ''//achar(27)//'[0m'
write(0,*) ''//achar(27)//'[0;38;5;033m' , 'a2(j) = ' , (a2(j),j=1,3) , ''//achar(27)//'[0m'
write(0,*) ''//achar(27)//'[0;38;5;033m' , 'a3(j) = ' , (a3(j),j=1,3) , ''//achar(27)//'[0m'
write(0,*) ' '










! interperet k-point index argument
nwk=nint(xnwk)
if (kpoint.gt.nwk) then
   write(0,*) '*** error - selected k=',kpoint,' > max k=',nwk
   stop
endif

! interperet band index argument
nband=nint(xnband)
if (band.gt.nband) then
   write(0,*) '*** error - selected band=',band,' > max band=',nband
   stop
endif

! allocate memory for "cener == band en
allocate(occ(nband))
allocate(cener(nband))










!!$*   POST PROCESSING FROM LATTICE AND ENCUT DATA

!!$* compute reciprocal lattice
call vcross(a2xa3,a2,a3)
Vcell=dot_product(a1,a2xa3)
a3mag=dsqrt(dot_product(a3,a3))
call vcross(b1,a2,a3)
call vcross(b2,a3,a1)
call vcross(b3,a1,a2)
!!$* normalize reciprocal lattice
b1=2.*pi*b1/Vcell
b2=2.*pi*b2/Vcell
b3=2.*pi*b3/Vcell
!!$* write out reciprocal lattice to 'GCOEFF_MOD.txt'
write(11,*) (sngl(b1(j)),j=1,3)
write(11,*) (sngl(b2(j)),j=1,3)
write(11,*) (sngl(b3(j)),j=1,3)
!! get lattice magnitudes
b1mag=dsqrt(b1(1)**2+b1(2)**2+b1(3)**2)
b2mag=dsqrt(b2(1)**2+b2(2)**2+b2(3)**2)
b3mag=dsqrt(b3(1)**2+b3(2)**2+b3(3)**2)

!!$ ---------- PRINT STATEMENTS -------------
write(0,*) ''//achar(27)//'[0;38;5;054m' , '[2] --- Reciprocals', ''//achar(27)//'[0m'
write(0,*) ''//achar(27)//'[0;38;5;162m' , 'b1(j) = ',(b1(j),j=1,3) , ''//achar(27)//'[0m'
write(0,*) ''//achar(27)//'[0;38;5;162m' , 'b2(j) = ',(b2(j),j=1,3) , ''//achar(27)//'[0m'
write(0,*) ''//achar(27)//'[0;38;5;162m' , 'b3(j) = ',(b3(j),j=1,3) , ''//achar(27)//'[0m'
write(0,*) ''//achar(27)//'[0;38;5;162m' , 'b1mag = ',b1mag , ''//achar(27)//'[0m'
write(0,*) ''//achar(27)//'[0;38;5;162m' , 'b2mag = ',b2mag , ''//achar(27)//'[0m'
write(0,*) ''//achar(27)//'[0;38;5;162m' , 'b3mag = ',b3mag , ''//achar(27)//'[0m'
write(0,*) ' '

phi12=acos((b1(1)*b2(1)+b1(2)*b2(2)+b1(3)*b2(3))/(b1mag*b2mag)) !! B1-B2 angle
call vcross(vtmp,b1,b2) !! a3*
vmag=dsqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2) !! norm( A3* )
sinphi123=(b3(1)*vtmp(1)+b3(2)*vtmp(2)+b3(3)*vtmp(3))/(vmag*b3mag) !! cos( A3*-B3 )
nb1maxA=(dsqrt(ecut*c)/(b1mag*abs(sin(phi12))))+1 !! 1 + K0 / norm(b1) / ||sin(B1-B2)||
nb2maxA=(dsqrt(ecut*c)/(b2mag*abs(sin(phi12))))+1 !! 1 + K0 / norm(b2) / ||sin(B1-B2)||
nb3maxA=(dsqrt(ecut*c)/(b3mag*abs(sinphi123)))+1  !! 1 + K0 / norm(b1) / ||cos(A3*-B3)||
npmaxA=nint(4.*pi*nb1maxA*nb2maxA*nb3maxA/3.) !! 4pi/3 * nb1 * nb2 * nb3

!!$ ---------- PRINT STATEMENTS -------------
write(0,*) ''//achar(27)//'[0;38;5;054m' , '[3] --- PostProcess', ''//achar(27)//'[0m'
write(0,*) ''//achar(27)//'[0;38;5;162m' , 'phi12 = ',phi12 , ''//achar(27)//'[0m'
write(0,*) ''//achar(27)//'[0;38;5;162m' , 'sinphi123 = ',sinphi123 , ''//achar(27)//'[0m'
write(0,*) ''//achar(27)//'[0;38;5;162m' , 'nb1maxA = ',nb1maxA , ''//achar(27)//'[0m'
write(0,*) ''//achar(27)//'[0;38;5;162m' , 'nb2maxA = ',nb2maxA , ''//achar(27)//'[0m'
write(0,*) ''//achar(27)//'[0;38;5;162m' , 'nb3maxA = ',nb3maxA , ''//achar(27)//'[0m'
write(0,*) ''//achar(27)//'[0;38;5;162m' , 'npmaxA = ',npmaxA , ''//achar(27)//'[0m'
write(0,*) ''//achar(27)//'[0;38;5;162m' , '-------', ''//achar(27)//'[0m'

phi13=acos((b1(1)*b3(1)+b1(2)*b3(2)+b1(3)*b3(3))/(b1mag*b3mag))
call vcross(vtmp,b1,b3)
vmag=dsqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
sinphi123=(b2(1)*vtmp(1)+b2(2)*vtmp(2)+b2(3)*vtmp(3))/(vmag*b2mag)
phi123=abs(asin(sinphi123))
nb1maxB=(dsqrt(ecut*c)/(b1mag*abs(sin(phi13))))+1
nb2maxB=(dsqrt(ecut*c)/(b2mag*abs(sinphi123)))+1
nb3maxB=(dsqrt(ecut*c)/(b3mag*abs(sin(phi13))))+1
npmaxB=nint(4.*pi*nb1maxB*nb2maxB*nb3maxB/3.)

!!$ ---------- PRINT STATEMENTS -------------
write(0,*) ''//achar(27)//'[0;38;5;162m' , 'phi13 = ',phi13 , ''//achar(27)//'[0m'
write(0,*) ''//achar(27)//'[0;38;5;162m' , 'sinphi123 = ',sinphi123 , ''//achar(27)//'[0m'
write(0,*) ''//achar(27)//'[0;38;5;162m' , 'nb1maxB = ',nb1maxB , ''//achar(27)//'[0m'
write(0,*) ''//achar(27)//'[0;38;5;162m' , 'nb2maxB = ',nb2maxB , ''//achar(27)//'[0m'
write(0,*) ''//achar(27)//'[0;38;5;162m' , 'nb3maxB = ',nb3maxB , ''//achar(27)//'[0m'
write(0,*) ''//achar(27)//'[0;38;5;162m' , 'npmaxB = ',npmaxB , ''//achar(27)//'[0m'
write(0,*) ''//achar(27)//'[0;38;5;162m' , '-------', ''//achar(27)//'[0m'

phi23=acos((b2(1)*b3(1)+b2(2)*b3(2)+b2(3)*b3(3))/(b2mag*b3mag))
call vcross(vtmp,b2,b3)
vmag=dsqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
sinphi123=(b1(1)*vtmp(1)+b1(2)*vtmp(2)+b1(3)*vtmp(3))/(vmag*b1mag)
phi123=abs(asin(sinphi123))
nb1maxC=(dsqrt(ecut*c)/(b1mag*abs(sinphi123)))+1
nb2maxC=(dsqrt(ecut*c)/(b2mag*abs(sin(phi23))))+1
nb3maxC=(dsqrt(ecut*c)/(b3mag*abs(sin(phi23))))+1 
npmaxC=nint(4.*pi*nb1maxC*nb2maxC*nb3maxC/3.)

!!$ ---------- PRINT STATEMENTS -------------
write(0,*) ''//achar(27)//'[0;38;5;162m' , 'phi23 = ',phi23 , ''//achar(27)//'[0m'
write(0,*) ''//achar(27)//'[0;38;5;162m' , 'sinphi123 = ',sinphi123 , ''//achar(27)//'[0m'
write(0,*) ''//achar(27)//'[0;38;5;162m' , 'nb1maxC = ',nb1maxC , ''//achar(27)//'[0m'
write(0,*) ''//achar(27)//'[0;38;5;162m' , 'nb2maxC = ',nb2maxC , ''//achar(27)//'[0m'
write(0,*) ''//achar(27)//'[0;38;5;162m' , 'nb3maxC = ',nb3maxC , ''//achar(27)//'[0m'
write(0,*) ''//achar(27)//'[0;38;5;162m' , 'npmaxC = ',npmaxC , ''//achar(27)//'[0m'
write(0,*) ''//achar(27)//'[0;38;5;162m' , '-------', ''//achar(27)//'[0m'

!!$* identify number of waves along each r-lattice vector
nb1max=max0(nb1maxA,nb1maxB,nb1maxC)
nb2max=max0(nb2maxA,nb2maxB,nb2maxC)
nb3max=max0(nb3maxA,nb3maxB,nb3maxC)

!! 2x to handle two component spinors
npmax=2*min0(npmaxA,npmaxB,npmaxC)

!!$*   Find desired wavefunction
irec=3+(kpoint-1)*(nband+1)

!!$ ---------- PRINT STATEMENTS -------------
write(0,*) ''//achar(27)//'[0;38;5;162m' , 'nb1max = ',nb1max , ''//achar(27)//'[0m'
write(0,*) ''//achar(27)//'[0;38;5;162m' , 'nb2max = ',nb2max , ''//achar(27)//'[0m'
write(0,*) ''//achar(27)//'[0;38;5;162m' , 'nb3max = ',nb3max , ''//achar(27)//'[0m'
write(0,*) ''//achar(27)//'[0;38;5;162m' , 'npmax = ',npmax , ''//achar(27)//'[0m'
write(0,*) ''//achar(27)//'[0;38;5;162m' , 'irec = ',irec , ''//achar(27)//'[0m'
write(0,*) ' '










!!$* large allocation
allocate (igall(3,npmax))
allocate (coeff(npmax))










!!$*   Find desired wavefunction
read(unit=10,rec=irec) xnplane,(wk(i),i=1,3), (cener(iband),occ(iband),iband=1,nband)
nplane=nint(xnplane)

!!$ ---------- PRINT STATEMENTS -------------
write(0,*) ''//achar(27)//'[0;38;5;021m' , '[4] --- Header Variables' , ''//achar(27)//'[0m'
write(0,*) ''//achar(27)//'[0;38;5;021m' , '    read(unit=10,rec=irec) xnplane,(wk(i),i=1,3),' , &
'(cener(iband),occ(iband),iband=1,nband)' , ''//achar(27)//'[0m'
write(0,*) ''//achar(27)//'[0;38;5;033m' , 'xnplane = ' , xnplane , ''//achar(27)//'[0m'
write(0,*) ''//achar(27)//'[0;38;5;033m' , 'wk = ' , (wk(i),i=1,3) , ''//achar(27)//'[0m'
write(0,*) ''//achar(27)//'[0;38;5;033m' , 'cener(1:nband) = ',(cener(iband),iband=1,nband), ''//achar(27)//'[0m'
write(0,*) ''//achar(27)//'[0;38;5;033m' , 'occ(1:nband) = ',(occ(iband),iband=1,nband), ''//achar(27)//'[0m'

!!$ ---------- "OLD" PRINT STATEMENTS -------------
!do idx=1,nband/2
!    write(0,*) ''//achar(27)//'[0;38;5;027m' ,(cener(iband),iband=2*idx-1,2*idx), ''//achar(27)//'[0m'
!end do
!write(0,*) ''//achar(27)//'[0;38;5;033m' , 'occ(1:nband) = ', ''//achar(27)//'[0m'
!do idx=1,nband/2
!    write(0,*) ''//achar(27)//'[0;38;5;027m' ,(occ(iband),iband=2*idx-1,2*idx), ''//achar(27)//'[0m'
!end do





!!$* write more details to 'GCOEFF_MOD.txt' header
write(11,*) kpoint,band
write(11,*) (sngl(wk(j)),j=1,3)
write(11,*) cener(band),occ(band)






!!$*   Calculate plane waves
ncnt=0
do ig3=0,2*nb3max
   ig3p=ig3
   if (ig3.gt.nb3max) ig3p=ig3-2*nb3max-1
   do ig2=0,2*nb2max
      ig2p=ig2
      if (ig2.gt.nb2max) ig2p=ig2-2*nb2max-1
      do ig1=0,2*nb1max
         ig1p=ig1
         if (ig1.gt.nb1max) ig1p=ig1-2*nb1max-1
         do j=1,3
            sumkg(j) = &
              (wk(1)+ig1p)*b1(j) + &
              (wk(2)+ig2p)*b2(j) + &
              (wk(3)+ig3p)*b3(j)
         enddo
         gtot=sqrt(dot_product(sumkg,sumkg))
         etot=gtot**2/c
         if (etot.lt.ecut) then
          ncnt=ncnt+1
          igall(1,ncnt)=ig1p
          igall(2,ncnt)=ig2p
          igall(3,ncnt)=ig3p
         else
          write(0,*) etot
         end if
      enddo
   enddo
enddo







! Acknowledge possible spinor failure
if (2*ncnt.ne.nplane) then
   write(0,*) '*** error - computed 2*ncnt=',2*ncnt, &
        ' != input nplane=',nplane
   stop
endif





irec=irec+band
read(unit=10,rec=irec) (coeff(iplane), iplane=1,nplane)
!!$*   output G value and coefficients
do iplane=1,nplane
   write(11,570) (igall(j,mod(iplane-1,ncnt)+1),j=1,3), coeff(iplane)
570 format(3i6,'  ( ',g14.6,' , ',g14.6,' )')     
enddo










! --------------- real-space calculations ---------------
do iz=0,2*nb3max
   z=dble(iz)/dble(1+2*nb3max)
   xyz(3)=z
   csum1=cmplx(0.,0.)
   csum2=cmplx(0.,0.)
   do iplane=1,ncnt
      ig=igall(:,iplane)
      wkpg=wk+ig
      csum1=csum1+coeff(iplane)* &
           cdexp(2.*pi*cmplx(0.,1.)*dot_product(wkpg,xyz))
      csum2=csum2+coeff(ncnt+iplane)* &
           cdexp(2.*pi*cmplx(0.,1.)*dot_product(wkpg,xyz))
   enddo
   csum1=csum1/dsqrt(Vcell)
   csum2=csum2/dsqrt(Vcell)
!!$ output z*a3mag for units of Angstroms
   write(6,*) sngl(z),sngl(real(csum1)),sngl(dimag(csum1)),&
        sngl(real(csum2)),sngl(dimag(csum2))
enddo




!! closing out...
deallocate(igall)
deallocate(coeff)
stop
end program




























!!$* -----------------------------------------------------------------------
!!$* -------------------------- sub-routines -------------------------------
!!$* -----------------------------------------------------------------------







!!$* compute vector cross-product
subroutine vcross(a,b,c)
  implicit real*8(a-h,o-z)
  dimension a(3),b(3),c(3)
  
  a(1)=b(2)*c(3)-b(3)*c(2)
  a(2)=b(3)*c(1)-b(1)*c(3)
  a(3)=b(1)*c(2)-b(2)*c(1)
  return
end subroutine vcross








!!$* parse command line arguments
subroutine parse(filename,kpoint,band,x,y)
character*75 filename
integer band,kpoint
real*8 x,y
character*20 option,value
integer iarg,narg,ia
iarg=iargc()
nargs=iarg/2
filename="WAVECAR"
kpoint = 1
band = 1
x = 0.
y = 0.
if(iarg.ne.2*nargs) then   ! must have an even number of inputs & flags
   call help
endif
do ia=1,nargs
   call getarg(2*ia-1,option)    ! idx1 = input index  ,  idx2 = the idx1 input element
   call getarg(2*ia,value)
   if(option == "-f") then
      read(value,*) filename     ! 'value' serves as the unit
   else if(option == "-k") then
      read(value,*) kpoint
   else if(option == "-b") then
      read(value,*) band
   else if(option == "-x") then
      read(value,*) x
   else if(option == "-y") then
      read(value,*) y
   else if(option =="-h") then
      call help
   else
      call help
   endif
enddo
return
end subroutine parse








!!$* UH OH! This prints a help statment for the user
subroutine help
  write(6,*) 'syntax: WaveTransSpinor -f file -k k-point -b band -x coord -y coord'
  write(6,*) 'defaults: -f WAVECAR -k 1 -b 1 -x 0.0 -y 0.0'
  write(6,*) 'inputs: x and y are direct coordinates on axes a1 and a2'
  write(6,*) 'output: spinor phi(x,y,z), chi(x,y,z) with z direct coordinate on a3 axis'
  stop
end subroutine help


