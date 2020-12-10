
MODULE vtylib__POSCAR

  ! --- 0 --- workspace ---------------------------------------------- !
  IMPLICIT none
  REAL*8, PARAMETER :: pi=4.*atan(1.)
  REAL*8 :: c
  data c/0.262465831d0/



  ! --- POSCAR data type
  TYPE poscarFile
	CHARACTER*72 :: fileName , header
    INTEGER :: ionNO , unitNO , lines
    REAL*8 :: volume
    REAL*8,DIMENSION(3) :: a1 , a2 , a3
	REAL*8,DIMENSION(3) :: b1 , b2 , b3
	REAL*8,ALLOCATABLE :: ionD(:,:) , ionC(:,:)
	CHARACTER*2,ALLOCATABLE :: species
  end TYPE poscarFile





  CONTAINS




! --- pull real values from line ------------------------------------- !
SUBROUTINE pullThreeReals( strIN , vecOUT )
  ! -- 0 - workspace
  IMPLICIT none
  CHARACTER*72,INTENT(in) :: strIN
  REAL*8,DIMENSION(3),INTENT(out) :: vecOUT
  CHARACTER*72 :: str1TMP , str2TMP
  INTEGER :: idx1 , idx2
  ! get indices
  str1TMP = trim(adjustl(strIN))
  idx1 = index(str1TMP," ")
  idx2 = index(str1TMP(1:len_trim(str1TMP))," ",.TRUE.)
  ! acquire vecOUT(1)
  str2TMP = trim(adjustl(str1TMP(1:idx1)))
  read(str2TMP,'(F19.15)') vecOUT(1)
  ! acquire vecOUT(2)
  str2TMP = trim(adjustl(str1TMP(idx1:idx2)))
  read(str2TMP,"(F19.15)") vecOUT(2)
  ! acquire vecOUT(3)
  str2TMP = trim(adjustl(str1TMP(idx2:72)))
  read(str2TMP,"(F19.15)") vecOUT(3)
  ! done ..
  return
end SUBROUTINE pullThreeReals







! --- vector cross-product ------------------------------------------- !
SUBROUTINE vcrossPOS( a , b , c )
    ! --- 0 --- workspace
    IMPLICIT none
    REAL*8,DIMENSION(3),INTENT(in) :: a , b
    REAL*8,DIMENSION(3),INTENT(out) :: c
    ! --- 1 --- code
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
    return
end SUBROUTINE vcrossPOS




SUBROUTINE pullPosSpecies( inSTR , posFile )
    ! --- 0 --- workspace
    IMPLICIT none
    CHARACTER*72,INTENT(inout) :: inSTR
    TYPE(poscarFile),INTENT(incd ouout) :: arraySTR
    INTEGER :: count=0 , j
    CHARACTER*72 :: tmpSTR=""
    LOGICAL :: blankSTATE=.TRUE.
    inSTR = " "//inSTR(1:71)
    do j = 2,len_trim(inSTR)+1
    	! increase tic
    	if((blankSTATE).and.(inSTR(j:j).ne." ")) then
    		count = count + 1
    		blankSTATE = .FALSE.
    	else if((.not.blankSTATE).and.(inSTR(j:j).eq." ")) then
    		tmpSTR = trim(tmpSTR)//" "//adjustl(inSTR(j-2:j-1))
    		blankSTATE = .TRUE.
    	endif
    	! end
    	write(0,*) inSTR(j:j) , " - " , blankSTATE
    enddo
    write(0,*) "count: " , count
    write(0,*) "tmpSTR: __" , trim(adjustl(tmpSTR)) , "__"
end SUBROUTINE pullPosSpecies






! --- read poscar ---------------------------------------------------- !
SUBROUTINE readPOSCAR( posFile )
    ! workspace
  	IMPLICIT none
  	TYPE(poscarFile),INTENT(inout) :: posFile
  	INTEGER :: iost , i , j
  	CHARACTER*72 :: testSTR , test1 , test2 , test3
  	REAL*8 :: scaling
  	REAL*8,DIMENSION(3) :: vecTMP
	! set initial values
	if (posFile%unitNO.lt.9) then
		write(0,*) "** input poscar file unit was below 9 set to 9."
		posFile%unitNO = 9
	else if (posFile%unitNo.gt.99) then
		write(0,*) "** input poscar file unit was above 99 set to 99."
		posFile%unitNO = 99
	endif
	! filename adjustment
	posFile%fileName = trim( posFile%fileName )
	posFile%volume = 0.
	! check the file exists
	open( unit=posFile%unitNO , &
        file=posFile%fileName , &
        iostat=iost , &
 		status='old' )
	if (iost.ne.0) then
		write(6,*) 'error opening file'
		write(0,*) 'Fail opening POSCAR named: ' , posFile%fileName
		write(0,*) 'Error code: ' , iost
		write(0,*) 'https://www.ibm.com/support/&
		knowledgecenter/en/SS3KZ4_9.0.0/com.ibm.xlf111.bg.doc/&
		xlflr/iostatvalues.htm'
		close(unit=posFile%unitNO)
		stop
 	endif
 	! 1st-pass reading
 	do j = 1,10000
		read(posFile%unitNO,'(A)',END=900) testSTR
		testSTR = trim(adjustl(testSTR))
		if (j.eq.1) posFile%header = testSTR
		if (j.eq.2) read(testSTR,'(F19.15)') scaling
		if (j.eq.3) CALL pullThreeReals( testSTR , posFile%a1 )
		if (j.eq.4) CALL pullThreeReals( testSTR , posFile%a2 )
		if (j.eq.5) CALL pullThreeReals( testSTR , posFile%a3 )
		if (j.eq.6) CALL pullPosSpecies( testSTR , posFile%species )
		write(0,*) j , " " , testSTR
	enddo
	900 posFile%lines = j-1
	! reset file
	close( unit=posFile%unitNO )
	open( unit=posFile%unitNO , &
        file=posFile%fileName , &
        iostat=iost , &
 		status='old' )
	! 2nd-pas reading
	do j = 1,posFile%lines
		read(posFile%unitNO,'(A)') testSTR
	enddo
	! adjust lattice by scaling, and compute volume
	posFile%a1 = scaling * posFile%a1
	posFile%a2 = scaling * posFile%a2
	posFile%a3 = scaling * posFile%a3
	! acquire volume
	CALL vcrossPOS( posFile%a2 , posFile%a3 , vecTMP )
	posFile%volume = dot_product( posFile%a1 , vecTMP )
	! compute reciprocals
	CALL vcrossPOS( posFile%a2 , posFile%a3 , posFile%b1 )
	CALL vcrossPOS( posFile%a3 , posFile%a1 , posFile%b2 )
	CALL vcrossPOS( posFile%a1 , posFile%a2 , posFile%b3 )
	posFile%b1 = (2.*pi/posFile%volume) * posFile%b1
	posFile%b2 = (2.*pi/posFile%volume) * posFile%b2
	posFile%b3 = (2.*pi/posFile%volume) * posFile%b3
	write(0,*) " "
	! done ...
 	return
end SUBROUTINE readPOSCAR





! --- print user data ------------------------------------------------ !
SUBROUTINE printPOSCAR( posFile )
  ! workspace
  IMPLICIT none
  TYPE(poscarFile),INTENT(inout) :: posFile
  CHARACTER*72 :: title , bold , reg , reset
  ! formatStrings
  title = ""//achar(27)//"[1;3;38;5;057m"
  bold = ""//achar(27)//"[1;3;38;5;141m"
  reg = ""//achar(27)//"[0;3;38;5;141m"
  reset = ""//achar(27)//"[0m"
  ! print user data
  write(0,"(A)") trim(title)//"----- POSCAR -----"
  write(0,"(A)") &
    trim(bold)//"fileName: "//&
    trim(reg)//trim(posFile%fileName)
  write(0,"(A,3F16.9)") trim(bold)//&
    "a1 "//trim(reg) , posFile%a1
  write(0,"(A,3F16.9)") trim(bold)//&
    "a2 "//trim(reg) , posFile%a2
  write(0,"(A,3F16.9)") trim(bold)//&
    "a3 "//trim(reg) , posFile%a3
  write(0,"(A,F16.9,A)") trim(bold)//&
    "volume: "//trim(reg) , posFile%volume , " A^3"
  write(0,"(A,3F16.9)") trim(bold)//&
    "b1 "//trim(reg) , posFile%b1
  write(0,"(A,3F16.9)") trim(bold)//&
    "b2 "//trim(reg) , posFile%b2
  write(0,"(A,3F16.9)") trim(bold)//&
    "b3 "//trim(reg) , posFile%b3
  write(0,*) " "
  ! done ...
  return
end SUBROUTINE printPOSCAR







end MODULE vtylib__POSCAR
