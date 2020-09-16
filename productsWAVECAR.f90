! -------1---------2---------3---------4---------5---------6---------7 !
! -------------------------------------------------------------------- !
! Hello!
! This bastardized script also originated from 'wavetransSpinor.f90'
! code. It computes Bloch Wave overlaps using a real-space integration.
! You can specify the central Bloch Wave using -k and -b flags.
! Further you can choose the range of bands/kpoints you evaluate -B -K.
! (This is important, since these integrals can get expensive...)
!
! For example, (after you've compiled 'productsWAVECAR.exe'):
!    ./parseWAVECAR -f WAVECAR -k 80 -b 45 -K 70-80 -B 40-60
!
! (Defaults are -f WAVECAR -k 1 -b 1 -K 1-5 -B 1-5)
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
!  September 16th, 2020
! -------------------------------------------------------------------- !


! --- 0 --- workspace
! file units:   10-WAVECAR  17-outputFile
! recall that 0,5,6 are reserved by Fortran
implicit real*8 (a-h, o-z)
complex*8 csum1 , csum2 , overlapValue
dimension a1(3) , a2(3) , a3(3)
dimension b1(3) , b2(3) , b3(3)
dimension a2xa3(3) , sumkg(3) , vtmp(3)
dimension wk(3) , wkC(3) , xyz(3) , xyzC(3) , wkpg(3) , ig(3)
integer kpoint ,Kmin , Kmax , band , Bmin , Bmax , totalRealTerms
character*75 fileName
character*75 samplingFileName
character*75 outputFileName
character*4 nbandString
dimension wk0(3)
complex*8, allocatable :: coeff(:) , coeff0(:)
complex*8, allocatable :: Rcoeff(:,:) , Rcoeff0(:,:)
integer, allocatable :: igall(:,:) , igall0(:,:)
real*8, allocatable :: sampling(:,:)
logical writeDirect , quietComments , sFlagged
! constant   c = 2m/hbar**2  (in 1/eV Ang^2, altered to match VASP)
data c/0.262465831d0/
pi=4.*atan(1.)










! --- 1 --- parse arguments
call parse(fileName,kpoint,Kmin,Kmax,band,Bmin,Bmax)










! --- 2 --- run WAVECAR file checks, and acquire header data
nrecl=24
open(unit=10,file=fileName,access='direct',recl=nrecl, &
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
open(unit=10,file=fileName,access='direct',recl=nrecl, &
     iostat=iost,status='old')
if (iost.ne.0) write(6,*) 'open error - iostat =',iost
! pull vasp variables (kpts,bands,encut,lattice)
read(unit=10,rec=2) xnwk , xnband , encut , &
(a1(j),j=1,3) , (a2(j),j=1,3) , (a3(j),j=1,3)
nwk=nint(xnwk)
nband=nint(xnband)
! adjust k,b selection ranges to new data
if(Kmax.gt.nwk) Kmax = nwk
if(Bmax.gt.nband) Bmax = nband
if(Kmin.gt.Kmax) Kmin = 1
if(Bmin.gt.Bmax) Bmin = 1
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
allocate ( coeff0(npmax) )
allocate ( igall(3,npmax) )
allocate ( igall0(3,npmax) )






! --- 4 --- set real-space sampling (and handle RPOINTS)
! sample frequencies
iXsample = 1*nb1max/2
iYsample = 1*nb2max/2
iZsample = 1*nb3max/8
iSampleCount=(iXsample+1)*(iYsample+1)*(iZsample+1)
! allocate memory
allocate ( sampling(iSampleCount,3) )
allocate ( Rcoeff(iSampleCount,2) )
allocate ( Rcoeff0(iSampleCount,2) )
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
write(0,*) "i = ",i
write(0,*) "iSampleCount = ",iSampleCount










! --- 5 --- acquire seed wavefunction
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
! CONVERT seed to REALSPACE (can be expensive)
write(0,*) "Calculating real-space seed wavefunction"
call convertToRealSpace( Rcoeff , coeff , npmax , iSampleCount , &
  sampling , xyz , nplane , wk , igall , Vcell )
! compute norm
waveNorm = sqrt(   dot_product( Rcoeff(:,1) , Rcoeff(:,1) ) + &
dot_product( Rcoeff(:,2) , Rcoeff(:,2) )   )
write(0,*) "waveNorm = ",waveNorm






! --- 6 --- write overlaps
! create file name
outputFileName=""//trim(fileName)//"_overlaps"
call makeFileName( outputFileName , kpoint , nwk , band , nband )
write(outputFileName,'(A)') ''//trim(outputFileName)//'.dat'
open(unit=17,file=outputFileName)
! write header
write(17,'(A)') 'Wavefunction overlap matrix elements:   < u_mq | u_nk >'
write(17,'(A)') ' '
write(17,'(A,I4.1,A,I4.1)') 'band:  n = ',band,' of ',nband
write(17,'(A,I4.1,A,I4.1)') 'k_idx: k = ',kpoint,' of ',nwk
write(17,'(A,3F14.9)') '  k_vec_direct    =',(wk(j),j=1,3)
wkC(1)=b1(1)*wk(1)+b2(1)*wk(2)+b3(1)*wk(3)
wkC(2)=b1(2)*wk(1)+b2(2)*wk(2)+b3(2)*wk(3)
wkC(3)=b1(3)*wk(1)+b2(3)*wk(2)+b3(3)*wk(3)
write(17,'(A,3F14.9,A)') '  k_vec_rec(A^-1) =',(wkC(j),j=1,3)
write(17,'(A)') ' '
write(17,'(A)') 'B =      (1/Angstrom)'
write(17,*) 'b1 = ' , (sngl(b1(j)),j=1,3)
write(17,*) 'b2 = ' , (sngl(b2(j)),j=1,3)
write(17,*) 'b3 = ' , (sngl(b3(j)),j=1,3)
write(17,'(A)') ' '
write(17,'(A,I4.1,A,I4.1)') 'Selected m range: ',Bmin,'   - ',Bmax
write(17,'(A,I4.1,A,I4.1)') 'Selected q range: ',Kmin,'   - ',Kmax
write(17,'(A)') ' '
write(17,'(A)') 'The four columns are real (x), imaginary, (y), &
                absolute value (r), and phase (theta;-pi,pi).'
write(17,'(A)') ' '
! BIG LOOPS BROTHER
do iq=Kmin,Kmax
  ! read new header data
  irec=3+(iq-1)*(nband+1)
  read(unit=10,rec=irec) xnplane,(wk0(i),i=1,3)
  nplane0=nint(xnplane)
  ! pull plane indices
  call yankPlanes( igall0 , npmax , &
  c,encut,nplane0,wk0,b1,b2,b3,nb1max,nb2max,nb3max )
  ! upper ribbon text
  write(17,*) " "
  write(17,'(A,I4.1)') "q = ",iq
  write(17,'(A,3F14.9)') "   vector q = " , ( wk0(j) , j=1,3 )
  wkC(1)=b1(1)*wk0(1)+b2(1)*wk0(2)+b3(1)*wk0(3)
  wkC(2)=b1(2)*wk0(1)+b2(2)*wk0(2)+b3(2)*wk0(3)
  wkC(3)=b1(3)*wk0(1)+b2(3)*wk0(2)+b3(3)*wk0(3)
  write(17,'(A,3F14.9)') "cartesian q = " , ( wkC(j) , j=1,3 )
  write(17,'(A,3F14.9)') "   vector q-k = " , ( wk0(j)-wk(j) , j=1,3 )
  wkC(1)=b1(1)*(wk0(1)-wk(1))+b2(1)*(wk0(2)-wk(2))+b3(1)*(wk0(3)-wk(3))
  wkC(2)=b1(2)*(wk0(1)-wk(1))+b2(2)*(wk0(2)-wk(2))+b3(2)*(wk0(3)-wk(3))
  wkC(3)=b1(3)*(wk0(1)-wk(1))+b2(3)*(wk0(2)-wk(2))+b3(3)*(wk0(3)-wk(3))
  write(17,'(A,3F14.9)') "cartesian q-k = " , ( wkC(j) , j=1,3 )
  ! LOOP bands
  do ib=Bmin,Bmax
    write(0,*) ""//achar(27)//"[0;1;38;5;064mstarting iteration ",&
              (ib-Bmin) + (Bmax-Bmin+1)*(iq-Kmin) + 1 , ' of ' ,&
              (Bmax-Bmin+1)*(Kmax-Kmin+1),""//achar(27)//"[0m"
    ! pull new coefficents
    irec=3+(iq-1)*(nband+1) + ib
    read(unit=10,rec=irec) (coeff0(iplane), iplane=1,nplane0)
    ! HEAVY convert to realspace
    write(0,*) "--> Calculating another real-space wavefunction"
    write(0,*) "---> iq =",iq ,"(Kmin =",Kmin , "  Kmax =",Kmax , ")"
    write(0,*) "---> ib =",ib ,"( Bmin =",Bmin , "  Bmax =",Bmax , ")"
    call convertToRealSpace( Rcoeff0 , coeff0 , npmax , iSampleCount &
                  , sampling , xyz , nplane0 , wk0 , igall0 , Vcell )
    ! compute norm
    waveNorm0 = sqrt(   dot_product( Rcoeff0(:,1) , Rcoeff0(:,1) ) &
      + dot_product( Rcoeff0(:,2) , Rcoeff0(:,2) )   )
    write(0,*)  "---> "//achar(27)//"[0;38;5;022m" ,&
                "norm0 = ",waveNorm0,""//achar(27)//"[0m"
    ! take inner product, and write result
    overlapValue = ( dot_product(   Rcoeff0(:,1) , Rcoeff(:,1) ) &
      + dot_product( Rcoeff0(:,2) , Rcoeff(:,2) )   ) &
      / waveNorm / waveNorm0
    tmpReal=REAL(REAL(overlapValue))
    tmpImag=REAL(AIMAG(overlapValue))
    write(17,570) "b = " , ib , overlapValue , &
                    abs(overlapValue) , atan2(tmpImag,tmpReal)
    570 format(a,i4.1,'  | ',g14.6,'  ',g14.6 ,' | ',f14.6,' ',f14.6)
  enddo
enddo

















! --- 8 --- exit program
! release memory resources
deallocate(igall)
deallocate(igall0)
deallocate(coeff)
deallocate(coeff0)
deallocate(Rcoeff)
deallocate(Rcoeff0)
deallocate(sampling)
! close files
close(10)
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
! workspace
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




! CONVERT to REALSPACE (can be expensive)
subroutine convertToRealSpace( Rcoeff , coeff , npmax , iSampleCount,&
sampling , xyz , nplane , wk , igall , Vcell )
! workspace
implicit real*8(a-h,o-z)
integer iSampleCount , nplane , npmax
complex*8, dimension(iSampleCount,2) :: Rcoeff
complex*8, dimension(npmax) :: coeff
integer, dimension(3,npmax) :: igall
integer, dimension(3) :: ig
real*8, dimension(iSampleCount,3) :: sampling
real*8, dimension(3) :: wk , xyz , wkpg
real*8 Vcell
complex*8 csum1 , csum2
! begin (potentially) massive loop
do i= 1,iSampleCount
    xyz(1) = sampling(i,1)
    xyz(2) = sampling(i,2)
    xyz(3) = sampling(i,3)
    csum1=cmplx(0.,0.)
    csum2=cmplx(0.,0.)
    do iplane=1,nplane/2
        ig=igall(:,iplane)
        wkpg = wk + ig
        csum1=csum1+coeff(iplane) * &
            cdexp(2.*pi*cmplx(0.,1.) * dot_product(wkpg,xyz))
        csum2=csum2+coeff(iplane+nplane/2) * &
            cdexp(2.*pi*cmplx(0.,1.) * dot_product(wkpg,xyz))
    enddo
    ! normalize output by volume
    csum1 = csum1/dsqrt(Vcell)
    csum2 = csum2/dsqrt(Vcell)
    ! write output
    Rcoeff(i,1) = csum1
    Rcoeff(i,2) = csum2
enddo
return
end subroutine convertToRealSpace













! parse command-line arguments
subroutine parse(filename,kpoint,Kmin,Kmax,band,Bmin,Bmax)
! workspace
integer kpoint,Kmin,Kmax,band,Bmin,Bmax
character*75 fileName , samplingFileName
logical writeDirect , quietComments , sFlagged
character*2, dimension(5) :: userFlags
character*20 option , value , val
integer iarg , narg , ia
! variable prep
iarg = iargc()
ia = 1
userFlags = (/"-b","-k","-f","-B","-K"/)
! user parameters
filename = "WAVECAR"
kpoint = 1
Kmin = 1
Kmax = 5
band = 1
Bmin = 1
Bmax = 5
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
  else if(option == "-K") then
    val = "0" // trim(adjustl(value)) // "0"
    if(index(val,"-").lt.2 .or. &
      index(val,"-").gt.len(val)) then
      call help
    else
      read(val(1:index(val,"-")-1),*) Kmin
      if(Kmin.eq.0) Kmin = 1
      read(val(index(val,"-")+1:len(val)),*) Kmax
      Kmax = Kmax/10
      if(Kmax.eq.0) Kmax = 10000
    endif
  else if(option == "-B") then
    val = "0" // trim(adjustl(value)) // "0"
    if(index(val,"-").lt.2 .or. &
      index(val,"-").gt.len(val)) then
      call help
    else
      read(val(1:index(val,"-")-1),*) Bmin
      if(Bmin.eq.0) Bmin = 1
      read(val(index(val,"-")+1:len(val)),*) Bmax
      Bmax = Bmax/10
      if(Bmax.eq.0) Bmax = 10000
    endif
  ! Unrecognized Flag!
  else
      call help
  endif
  ! increment counter
  ia=ia+1
enddo
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


