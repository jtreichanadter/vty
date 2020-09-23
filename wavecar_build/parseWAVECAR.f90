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
! (So for raw data, try 'tail -n +25 WAVECAR_k80_b45.dat' )
!
!
! OTHER FLAGS:
!
! -q ::: quiet; so no warnings or time indicators to outstream.
!
!
! -c ::: cartesian; writes coordinates for a wavefunction in cartesian
!                   coordinates. The default behavior (without the
!                   flag) is to write them in direct units.
!
! -s ::: sampling; looks for a file RPOINTS.
!                  If found (and with appropriate format), then the
!                  wavefunction will be evaluated at the points chosen
!                  in RPOINTS. Otherwise there is a default way I set
!                  the real-space sampling, and having the -s flag set
!                  in this case means an RPOINTS file will be created.
!                  (Also you can use another name, ex: -s MYPOINTS).
!
! RPOINTS format: The 1st line specifies DIRECT or CARTESIAN
!                 (defaults to direct, but if the first character is a
!                 'c' or 'C' then it will assume cartesian format).
!                 Every following line should be a triplet of floats.
!
! (Also, see compiling advice in the README.md file)
!
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


! --- 0 --- workspace
! file units:   10-WAVECAR   11-RPOINTS   17-outputFile
! recall that 0,5,6 are reserved by Fortran
implicit real*8 (a-h, o-z)
complex*8 csum1 , csum2
dimension a1(3) , a2(3) , a3(3)
dimension b1(3) , b2(3) , b3(3)
dimension a2xa3(3) , sumkg(3) , vtmp(3)
dimension wk(3) , wkC(3) , xyz(3) , xyzC(3) , wkpg(3) , ig(3)
integer kpoint , band , totalRealTerms
character*75 fileName
character*75 samplingFileName
character*75 outputFileName
character*4 nbandString
dimension wk0(3)
complex*8, allocatable :: coeff(:)
integer, allocatable :: igall(:,:)
real*8, allocatable :: sampling(:,:)
logical writeDirect , quietComments , sFlagged
! constant   c = 2m/hbar**2  (in 1/eV Ang^2, altered to match VASP)
data c/0.262465831d0/
pi=4.*atan(1.)







! --- 1 --- parse arguments
call parse( kpoint , band , fileName , samplingFileName , &
            writeDirect , quietComments , sFlagged )






! --- 2 --- run WAVECAR file checks, and acquire header data
nrecl=24
open(unit=10,file=filename,access='direct',recl=nrecl, &
     iostat=iost,status='old')
if (iost.ne.0) then
  write(6,*) 'open error - iostat =',iost
  write(0,*) 'Could not open WAVECAR file named ',fileName
  write(0,*) 'Error code ',iost
  write(0,*) 'check out https://www.ibm.com/support/&
  knowledgecenter/en/SS3KZ4_9.0.0/com.ibm.xlf111.bg.doc/&
  xlflr/iostatvalues.htm'
  close(unit=10)
  stop
endif
! pull header variables
read(unit=10,rec=1) xnrecl,xnspin,xnprec
close(unit=10)
nrecl=nint(xnrecl)
nspin=nint(xnspin)
nprec=nint(xnprec)
! only continue if correct precision identified
if(nprec.eq.45210) then
  write(0,*) '*** error - WAVECAR_double requires complex*16'
  write(0,*) '*** error - WAVECAR_double requires complex*16'
  stop
endif
if(nspin.eq.2) then
  write(0,*) '*** error - Not a spinor WAVECAR. ISPIN =',nspin
  write(0,*) 'For noncollinear magnetism, you actually need ISPIN=1'
  stop
endif
! re-open file
open(unit=10,file=filename,access='direct',recl=nrecl, &
     iostat=iost,status='old')
if (iost.ne.0) write(6,*) 'open error - iostat =',iost
! pull vasp variables (kpts,bands,encut,lattice)
read(unit=10,rec=2) xnwk , xnband , encut , &
(a1(j),j=1,3) , (a2(j),j=1,3) , (a3(j),j=1,3)
nwk=nint(xnwk)
nband=nint(xnband)
! more checks...
if (kpoint.gt.nwk) then
   write(0,*) '*** error - selected k=',kpoint,' > max k=',nwk
   stop
endif
if (band.gt.nband) then
   write(0,*) '*** error - selected band=',band,' > max band=',nband
   stop
endif






! --- 3 --- compute parsing information
! reciprocal lattice
call vcross(a2xa3,a2,a3)
Vcell=dot_product(a1,a2xa3)
a3mag=dsqrt(dot_product(a3,a3))
call vcross(b1,a2,a3)
call vcross(b2,a3,a1)
call vcross(b3,a1,a2)
b1=2.*pi*b1/Vcell
b2=2.*pi*b2/Vcell
b3=2.*pi*b3/Vcell
b1mag=dsqrt(b1(1)**2+b1(2)**2+b1(3)**2)
b2mag=dsqrt(b2(1)**2+b2(2)**2+b2(3)**2)
b3mag=dsqrt(b3(1)**2+b3(2)**2+b3(3)**2)
! A-vec cutoff
phi12=acos((b1(1)*b2(1)+b1(2)*b2(2)+b1(3)*b2(3))/(b1mag*b2mag))
call vcross(vtmp,b1,b2)
vmag=dsqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
sinphi123=(b3(1)*vtmp(1)+b3(2)*vtmp(2)+b3(3)*vtmp(3))/(vmag*b3mag)
nb1maxA=(dsqrt(encut*c)/(b1mag*abs(sin(phi12))))+1
nb2maxA=(dsqrt(encut*c)/(b2mag*abs(sin(phi12))))+1
nb3maxA=(dsqrt(encut*c)/(b3mag*abs(sinphi123)))+1
npmaxA=nint(4.*pi*nb1maxA*nb2maxA*nb3maxA/3.)
! B-vec cutoff
phi13=acos((b1(1)*b3(1)+b1(2)*b3(2)+b1(3)*b3(3))/(b1mag*b3mag))
call vcross(vtmp,b1,b3)
vmag=dsqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
sinphi123=(b2(1)*vtmp(1)+b2(2)*vtmp(2)+b2(3)*vtmp(3))/(vmag*b2mag)
phi123=abs(asin(sinphi123))
nb1maxB=(dsqrt(encut*c)/(b1mag*abs(sin(phi13))))+1
nb2maxB=(dsqrt(encut*c)/(b2mag*abs(sinphi123)))+1
nb3maxB=(dsqrt(encut*c)/(b3mag*abs(sin(phi13))))+1
npmaxB=nint(4.*pi*nb1maxB*nb2maxB*nb3maxB/3.)
! C-vec cutoff
phi23=acos((b2(1)*b3(1)+b2(2)*b3(2)+b2(3)*b3(3))/(b2mag*b3mag))
call vcross(vtmp,b2,b3)
vmag=dsqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
sinphi123=(b1(1)*vtmp(1)+b1(2)*vtmp(2)+b1(3)*vtmp(3))/(vmag*b1mag)
phi123=abs(asin(sinphi123))
nb1maxC=(dsqrt(encut*c)/(b1mag*abs(sinphi123)))+1
nb2maxC=(dsqrt(encut*c)/(b2mag*abs(sin(phi23))))+1
nb3maxC=(dsqrt(encut*c)/(b3mag*abs(sin(phi23))))+1
npmaxC=nint(4.*pi*nb1maxC*nb2maxC*nb3maxC/3.)
! formalize cutoffs
nb1max=max0(nb1maxA,nb1maxB,nb1maxC)
nb2max=max0(nb2maxA,nb2maxB,nb2maxC)
nb3max=max0(nb3maxA,nb3maxB,nb3maxC)
! x2 for handling two-component spinors
npmax=2*min0(npmaxA,npmaxB,npmaxC)
! allocate dynamic memory
allocate ( coeff(npmax) )
allocate ( igall(3,npmax) )








! --- 4 --- acquire seed wavefunction
! read header data
irec=3+(kpoint-1)*(nband+1)
read(unit=10,rec=irec) xnplane , (wk(i),i=1,3)
nplane=nint(xnplane)
! calculate planewave indices
call yankPlanes( igall , size(igall) , &
 c , encut , nplane , wk , b1 , b2 , b3 , nb1max , nb2max , nb3max )
! pull coefficents
irec=irec+band
read(unit=10,rec=irec) (coeff(iplane), iplane=1,nplane)




! --- 5 --- set real-space sampling (and handle RPOINTS)
! check RPOINTS
open(unit=11,file=samplingFileName,access='sequential', &
     iostat=iost,status='old')
close(11)

! CASE: read RPOINTS
if ((sFlagged).and.(iost==0)) then
  write(0,*) "Reading real-space sampling from ",samplingFileName
  ! open file
  open(unit=11,file=samplingFileName,access='sequential', &
    iostat=iost,status='old')
  ! read number of lines (to allocate memory)
  nlines = 0
  do
    read(10,*,iostat=iost)
    if (iost/=0) exit
    nlines = nlines + 1
  end do
  ! allocate memory
  allocate ( sampling(nlines-1,3) )
  ! STUB!!! [CHECK cartesian or direct]
! CASE: error, corrupt file??
else if ((sFlagged).and.(.not.(iost==6))) then
  write(6,*) "ERROR --- strange reading error of ",samplingFileName
  write(0,*) "iostat = ",iost
  stop
! CASE: create default unit-cell grid
else
  ! sample frequencies
  iXsample = 2*nb1max
  iYsample = 2*nb2max
  iZsample = 2*nb3max
  iSampleCount=(iXsample+1)*(iYsample+1)*(iZsample+1)
  ! allocate memory
  allocate ( sampling(iSampleCount,3) )
  ! fill in the sampling array
  i=0
  do iz=0,iZsample
  do iy=0,iYsample
  do ix=0,iXsample
    i=i+1
    sampling(i,1) = dble(ix) / dble(iXsample+1)
    sampling(i,2) = dble(iy) / dble(iYsample+1)
    sampling(i,3) = dble(iz) / dble(iZsample+1)
  enddo
  enddo
  enddo
  ! CASE: write to RPOINTS file
  if ((sFlagged).and.(iost==6)) then
    write(0,*) "Writing real-space sampling to ",samplingFileName
    ! STUB!!!
  endif
endif





! --- 6 --- compute real-space wavefunction
! create file name
outputFileName=fileName
call makeFileName( outputFileName , &
                kpoint , nwk , band , nband )
if (writeDirect) then
  write(outputFileName,'(A)') ''//trim(outputFileName)//'.dat'
else
  write(outputFileName,'(A)') ''//trim(outputFileName)//'_c.dat'
endif
! write header
open( unit=17 , file=outputFileName )
write(17,'(A)') 'Wavefunction coefficients real space'
write(17,'(A)') ' '
write(17,'(A)') 'Bloch wave details:'
write(17,'(A,I4.1,A,I4.1)') '  band:   n = ',band,' of ',nband
write(17,'(A,I4.1,A,I4.1)') '  k_idx:  k = ',kpoint,' of ',nwk
write(17,'(A,3F14.9)') '  k_vec_direct    =',(wk(j),j=1,3)
wkC(1)=b1(1)*wk(1)+b2(1)*wk(2)+b3(1)*wk(3)
wkC(2)=b1(2)*wk(1)+b2(2)*wk(2)+b3(2)*wk(3)
wkC(3)=b1(3)*wk(1)+b2(3)*wk(2)+b3(3)*wk(3)
write(17,'(A,3F14.9,A)') '  k_vec_rec(A^-1) =',(wkC(j),j=1,3)
write(17,'(A)') ' '
write(17,'(A,I8.1)') 'Total sampling points: ', iSampleCount
write(17,'(A,I6.1)') '  sampling x: ', iXsample
write(17,'(A,I6.1)') '  sampling y: ', iYsample
write(17,'(A,I6.1)') '  sampling z: ', iZsample
write(17,'(A)') ' '
write(17,'(A)') 'r =  n1 a1  +  n2 a2  +  n3 a3'
write(17,'(A)') 'A =      (Angstrom)'
write(17,*) 'a1 = ' , (sngl(a1(j)),j=1,3)
write(17,*) 'a2 = ' , (sngl(a2(j)),j=1,3)
write(17,*) 'a3 = ' , (sngl(a3(j)),j=1,3)
write(17,'(A)') ' '
write(17,'(A)') ' '
write(17,'(A)') ' '
write(17,'(A)') ' '
if (writeDirect) then
  write(17,'(A)') '      n1         n2         n3     &
  |       spin up (complex)         |      spin down (complex)       |'
else
  write(17,'(A)') '     x (A)      y (A)      z (A)   &
  |       spin up (complex)         |      spin down (complex)       |'
endif
write(11,'(A)') ' '
! - DO THE THING
do i= 1,iSampleCount
    xyz(1) = sampling(i,1)
    xyz(2) = sampling(i,2)
    xyz(3) = sampling(i,3)
    csum1=cmplx(0.,0.)
    csum2=cmplx(0.,0.)
    do iplane=1,nplane/2
        ig=igall(:,iplane)
        wkpg = wk + ig
        csum1=csum1+coeff(iplane)* &
            cdexp(2.*pi*cmplx(0.,1.)*dot_product(wkpg,xyz))
        csum2=csum2+coeff(iplane+nplane/2)* &
            cdexp(2.*pi*cmplx(0.,1.)*dot_product(wkpg,xyz))
    enddo
    ! normalize output by volume
    csum1 = csum1/dsqrt(Vcell)
    csum2 = csum2/dsqrt(Vcell)
    ! write output
    if (writeDirect) then
      write(17,571) (xyz(j),j=1,3) , csum1 , csum2
      571 format(3f11.5,'  &
      | ',g14.6,' ',g14.6,'   | ',g14.6,' ',g14.6,'  |')
    else
      xyzC(1)=a1(1)*xyz(1)+a2(1)*xyz(2)+a3(1)*xyz(3)
      xyzC(2)=a1(2)*xyz(1)+a2(2)*xyz(2)+a3(2)*xyz(3)
      xyzC(3)=a1(3)*xyz(1)+a2(3)*xyz(2)+a3(3)*xyz(3)
      write(17,572) (xyzC(j),j=1,3) , csum1 , csum2
      572 format(3f11.5,'  &
      | ',g14.6,' ',g14.6,'   | ',g14.6,' ',g14.6,'  |')
    endif
enddo















! --- 8 --- exit program
! release memory resources
deallocate(igall)
deallocate(coeff)
deallocate(sampling)
! close files
close(10)
close(11)
close(17)
! exit
stop
end program




































! -------------------------------------------------------------------- !
! -------------------------------------------------------------------- !
! --------------------------- subroutines ---------------------------- !
! -------------------------------------------------------------------- !
! -------------------------------------------------------------------- !




! Compute vector cross-product
subroutine vcross(a,b,c)
    implicit real*8(a-h,o-z)
    dimension a(3),b(3),c(3)
    a(1)=b(2)*c(3)-b(3)*c(2)
    a(2)=b(3)*c(1)-b(1)*c(3)
    a(3)=b(1)*c(2)-b(2)*c(1)
    return
end subroutine vcross








! Identify wavefunction plane-wave indices
subroutine yankPlanes(igall,nigall,c,encut,nplane,wk, &
b1,b2,b3,nb1max,nb2max,nb3max)
implicit real*8(a-h,o-z)
dimension wk(3),b1(3),b2(3),b3(3),sumkg(3)
dimension igall(3,nigall)
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
(wk(1)+ig1p)*b1(j) + (wk(2)+ig2p)*b2(j) + (wk(3)+ig3p)*b3(j)
         enddo
         gtot=sqrt(dot_product(sumkg,sumkg))
         etot=gtot**2/c
         if (etot.lt.encut) then
            ncnt=ncnt+1
            igall(1,ncnt)=ig1p
            igall(2,ncnt)=ig2p
            igall(3,ncnt)=ig3p
         end if
      enddo
   enddo
enddo
if (2*ncnt.ne.nplane) then
   write(0,*) '*** error - computed 2*ncnt=',2*ncnt, &
        ' != input nplane=',nplane
   stop
endif
end subroutine yankPlanes







! parse command-line arguments
subroutine parse( kpoint , band , filename , samplingFileName , &
                  writeDirect , quietComments , sFlagged )
! workspace
integer kpoint , band
character*75 fileName , samplingFileName
logical writeDirect , quietComments , sFlagged
character*2, dimension(6) :: userFlags
character*20 option , value , val
integer iarg , narg , ia
! variable prep
iarg = iargc()
ia = 1
userFlags = (/"-b","-k","-f","-s","-c","-q"/)
! user parameters
kpoint = 1
band = 1
filename = "WAVECAR"
samplingFileName = "RPOINTS"
writeDirect = .TRUE.
quietComments = .FALSE.
sFlagged = .FALSE.
! LOOP through arguments
do while ( ia < iarg )
  call getarg(ia,option)
  call getarg(ia+1,value)
  if (any( userFlags==value )) value = " "
  if (.not.(any( userFlags==value ))) ia = ia + 1
  ! kpoint
  if(option=="-k") then
    if (.not.(value==" ")) read(value,*) kpoint
  ! band
  else if(option=="-b") then
    if (.not.(value==" ")) read(value,*) band
  ! fileName
  else if(option=="-f") then
    if (.not.(value==" ")) read(value,*) fileName
  ! samplingFileName
  else if(option=="-s") then
    sFlagged = .TRUE.
    if (.not.(value==" ")) read(value,*) samplingFileName
  ! writeDirect
  else if(option=="-c") then
    writeDirect = .FALSE.
  ! quietComments
  else if(option=="-q") then
    quietComments = .TRUE.
  ! Unrecognized Flag!
  else
      call help
  endif
  ! increment counter
  ia=ia+1
enddo
! handle case if the final argument is missed
call getarg(iarg,option)
if (option=="-s") sFlagged=.TRUE.
if (option=="-c") writeDirect=.FALSE.
if (option=="-q") quietComments=.TRUE.
! done...
return
end subroutine parse







! print help message
subroutine help
  write(0,*) ''//achar(27)//"[1;3;38;5;125m"//"ERROR - Uh Oh!"
  write(0,*) ''//achar(27)//"[0;38;5;125m"//"syntax:"
  write(0,*)" parseWAVECAR.f90 arguments   -f file   -k kpoint   &
  -b band   -s samplingFile   -c use cartesian coordinates   &
  -q quiet (suppress commandline outstream)"
  write(0,*) "defaults:"
  write(0,*) " -f WAVECAR -k 1 -b 1 -f WAVECAR"
  write(0,*) "NOTE flags are case-sensitive"
  stop
end subroutine help








! Format wavefunction coefficent file name
subroutine makeFileName( outputFileName , &
                kpoint , nwk , band , nband )
integer band , kpoint
integer nband , nwk
character*75 outputFileName
if(nwk>999) then
    if(nband>999) then
        write(outputFileName,'(A,I4.4,A,I4.4)') &
''//trim(outputFileName)//'_k',kpoint,'_b',band
    else if(nband>99) then
        write(outputFileName,'(A,I4.4,A,I3.3)') &
''//trim(outputFileName)//'_k',kpoint,'_b',band
    else if(nband>9) then
        write(outputFileName,'(A,I4.4,A,I2.2)') &
''//trim(outputFileName)//'_k',kpoint,'_b',band
    else
        write(outputFileName,'(A,I4.4,A,I1.1)') &
''//trim(outputFileName)//'_k',kpoint,'_b',band
    endif
else if(nwk>99) then
    if(nband>999) then
        write(outputFileName,'(A,I3.3,A,I4.4)') &
''//trim(outputFileName)//'_k',kpoint,'_b',band
    else if(nband>99) then
        write(outputFileName,'(A,I3.3,A,I3.3)') &
''//trim(outputFileName)//'_k',kpoint,'_b',band
    else if(nband>9) then
        write(outputFileName,'(A,I3.3,A,I2.2)') &
''//trim(outputFileName)//'_k',kpoint,'_b',band
    else
        write(outputFileName,'(A,I3.3,A,I1.1)') &
''//trim(outputFileName)//'_k',kpoint,'_b',band
    endif
else if(nwk>9) then
    if(nband>999) then
        write(outputFileName,'(A,I2.2,A,I4.4)') &
''//trim(outputFileName)//'_k',kpoint,'_b',band
    else if(nband>99) then
        write(outputFileName,'(A,I2.2,A,I3.3)') &
''//trim(outputFileName)//'_k',kpoint,'_b',band
    else if(nband>9) then
        write(outputFileName,'(A,I2.2,A,I2.2)') &
''//trim(outputFileName)//'_k',kpoint,'_b',band
    else
        write(outputFileName,'(A,I2.2,A,I1.1)') &
''//trim(outputFileName)//'_k',kpoint,'_b',band
    endif
else
    if(nband>999) then
        write(outputFileName,'(A,I1.1,A,I4.4)') &
''//trim(outputFileName)//'_k',kpoint,'_b',band
    else if(nband>99) then
        write(outputFileName,'(A,I1.1,A,I3.3)') &
''//trim(outputFileName)//'_k',kpoint,'_b',band
    else if(nband>9) then
        write(outputFileName,'(A,I1.1,A,I2.2)') &
''//trim(outputFileName)//'_k',kpoint,'_b',band
    else
        write(outputFileName,'(A,I1.1,A,I1.1)') &
''//trim(outputFileName)//'_k',kpoint,'_b',band
    endif
endif
end subroutine makeFileName


