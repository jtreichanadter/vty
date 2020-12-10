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


MODULE vtylib__WAVECAR

  ! --- 0 --- workspace ---------------------------------------------- !
  IMPLICIT none
  REAL*8, PARAMETER :: pi=4.*atan(1.)
  REAL*8 :: c
  data c/0.262465831d0/



  ! --- 1 --- derived data types (data strucutres) ------------------- !


  ! --- User-proveded arguments
  TYPE userArgumentStruct
    CHARACTER*72 :: fileName_INPUT , fileName_RPOINTS
    INTEGER :: k , kMIN , kMAX
    INTEGER :: b , bMIN , bMAX
    LOGICAL :: readRPOINTS , writeRPOINTS
    LOGICAL :: isQuiet , isDirectRPOINTS
  end TYPE userArgumentStruct


  ! --- WAVECAR file metadata
  TYPE wavecarFile
    CHARACTER*72 :: fileName , magState
    INTEGER :: unitNo , kNO , bNO , enCut , npCAP
    INTEGER :: precision , spin , recordLength
    REAL*8 :: volume
    INTEGER,DIMENSION(3) :: nbMAX , npMAX
    REAL*8,DIMENSION(3) :: a1 , a2 , a3
    REAL*8,DIMENSION(3) :: b1 , b2 , b3
  end TYPE wavecarFile


  ! --- wave function class
  TYPE waveFunction
    TYPE(wavecarFile) :: src
    INTEGER :: kpoint , band , nplane
    REAL*8,DIMENSION(3) :: kvec
    COMPLEX*8,ALLOCATABLE :: uG(:,:) , uR(:,:)
    INTEGER,ALLOCATABLE :: idx(:,:)
  end TYPE waveFunction




CONTAINS






! --- vector cross-product ------------------------------------------- !
SUBROUTINE vcross( a , b , c )
    ! --- 0 --- workspace
    IMPLICIT none
    REAL*8,DIMENSION(3),INTENT(in) :: a , b
    REAL*8,DIMENSION(3),INTENT(out) :: c
    ! --- 1 --- code
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
    return
end SUBROUTINE vcross




! --- parse command-line arguments ----------------------------------- !
subroutine parse( user )
  ! workspace
  IMPLICIT none
  TYPE(userArgumentStruct),INTENT(inout) :: user
  CHARACTER*2,DIMENSION(8) :: userFlags
  CHARACTER*72 :: option , value , val
  INTEGER :: iarg , ia ,idx , tmp
  ! variable prep
  userFlags = (/"-b","-B","-k","-K","-f","-r","-c","-q"/)
  iarg = iargc()
  ia = 1
  ! user defaults
  user%b = 1
  user%bMIN = 1
  user%bMAX = 10000
  user%k = 1
  user%kMIN = 1
  user%kMAX = 10000
  user%fileName_INPUT = "WAVECAR"
  user%fileName_RPOINTS = "RPOINTS"
  user%isQuiet = .FALSE.
  user%isDirectRPOINTS = .TRUE.
  user%readRPOINTS = .FALSE.
  user%writeRPOINTS = .FALSE.
  ! LOOP through arguments
  do while ( ia .le. iarg )
    ! read commandline arguments
    call getarg(ia,option)
    call getarg(ia+1,value)
    ! clear value if it's a flag, else jump twice
    if ( .not.(any(userFlags==value)) ) ia=ia+1
    if ( any(userFlags==value) ) value = " "
    ! kpoint
    if(option=="-k") then
      if(value==" ") CALL userHelpMessage
      read(value,*) user%k
    ! kpoint range
    else if(option=="-K") then
      if(value==" ") CALL userHelpMessage
      val = "0" // trim(adjustl(value)) // "0"
      idx = index(val,"-")
      if(idx==0) CALL userHelpMessage
      read(val(1:idx-1),*) user%kMIN
      read(val(idx+1:len(val)),*) user%kMAX
      user%kMAX = user%kMAX/10
      if(user%kMIN.eq.0) user%kMIN = 1
      if(user%kMAX.eq.0) user%kMAX = 10000
      if(user%kMAX.le.user%kMIN) CALL userHelpMessage
    ! band
    else if(option=="-b") then
      if(value==" ") CALL userHelpMessage
      read(value,*) user%b
    ! band range
    else if(option=="-B") then
      if(value==" ") CALL userHelpMessage
      val = "0" // trim(adjustl(value)) // "0"
      idx = index(val,"-")
      if(idx==0) CALL userHelpMessage
      read(val(1:idx-1),*) user%bMIN
      read(val(idx+1:len(val)),*) user%bMAX
      user%bMAX = user%bMAX/10
      if(user%bMIN.eq.0) user%bMIN = 1
      if(user%bMAX.eq.0) user%bMAX = 10000
      if(user%bMAX.le.user%bMIN) CALL userHelpMessage
    ! fileName_INPUT
    else if(option=="-f") then
      if(value==" ") CALL userHelpMessage
      user%fileName_INPUT = value(1:index(value," "))
    ! fileName_RPOINTS
    else if(option=="-r") then
      user%readRPOINTS = .TRUE.
      user%fileName_RPOINTS = value(1:index(value," "))
    ! isDirectRPOINTS
    else if(option=="-c") then
      user%isDirectRPOINTS = .FALSE.
    ! quietComments
    else if(option=="-q") then
      user%isQuiet = .TRUE.
    ! Unrecognized Flag!
    else
        write(0,*) "Unidentified user flag!!"
        CALL userHelpMessage
    endif
    ! increment counter
    ia=ia+1
  enddo
  ! handle case if the final argument is missed
  call getarg( iarg , option )
  if (option=="-r") user%readRPOINTS = .TRUE.
  if (option=="-c") user%isDirectRPOINTS = .FALSE.
  if (option=="-q") user%isQuiet = .TRUE.
  ! done...
  return
end SUBROUTINE parse




! --- print user data ------------------------------------------------ !
SUBROUTINE printUser( userArgs )
  ! workspace
  IMPLICIT none
  TYPE(userArgumentStruct),INTENT(inout) :: userArgs
  CHARACTER*72 :: title , bold , reg , reset
  ! formatStrings
  title = ""//achar(27)//"[1;3;38;5;118m"
  bold = ""//achar(27)//"[1;3;38;5;192m"
  reg = ""//achar(27)//"[0;3;38;5;192m"
  reset = ""//achar(27)//"[0m"
  ! print user data
  write(0,"(A)") trim(title)//"----- USER INPUT -----"
  write(0,"(A)") &
    trim(bold)//"wavecar file name: "//&
    trim(reg)//trim(userArgs%fileName_INPUT)
  write(0,"(A,I5,A,I2,A,I5,A)") trim(bold)//&
    "band: "//trim(reg),&
    userArgs%b , &
    "   ( " , userArgs%bMIN , " - " , userArgs%bMAX , " )"
  write(0,"(A,I5,A,I2,A,I5,A)") trim(bold)//&
    "kpoint: "//trim(reg),&
    userArgs%k , &
    "   ( " , userArgs%kMIN , " - " , userArgs%kMAX , " )"
  write(0,"(A,A)") trim(bold)//&
    "RPOINTS fileName: "//trim(reg) , userArgs%fileName_RPOINTS
  write(0,"(A,L)") trim(bold)//&
    "run quietly: "//trim(reg) , userArgs%isQuiet
  write(0,"(A,L)") trim(bold)//&
    "RPOINTS will be read: "//trim(reg) , userArgs%readRPOINTS
  write(0,"(A,L)") trim(bold)//&
    "RPOINTS are direct: "//trim(reg) , userArgs%isDirectRPOINTS
  write(0,"(A,L,A)") trim(bold)//&
    "RPOINTS will be written: "//trim(reg) , userArgs%writeRPOINTS
  write(0,"(A)") reset
end SUBROUTINE printUser





! --- print user data ------------------------------------------------ !
SUBROUTINE printWAVECAR( waveFile )
  ! workspace
  IMPLICIT none
  TYPE(wavecarFile),INTENT(inout) :: waveFile
  CHARACTER*72 :: title , bold , reg , reset
  INTEGER :: j
  ! formatStrings
  title = ""//achar(27)//"[1;3;38;5;039m"
  bold = ""//achar(27)//"[1;3;38;5;153m"
  reg = ""//achar(27)//"[0;3;38;5;153m"
  reset = ""//achar(27)//"[0m"
  ! print wavecar data
  write(0,"(A)") trim(title)//"----- WAVECAR -----"
  write(0,"(A)") &
    trim(bold)//"file name: "//&
    trim(reg)//trim(waveFile%fileName)
  write(0,"(A)") &
    trim(bold)//"magnetic format: "//&
    trim(reg)//trim(waveFile%magState)
  write(0,"(A,I4,A)") &
    trim(bold)//"energy cutoff: "//&
    trim(reg) , waveFile%enCut , " eV"
  write(0,"(A,I4)") &
    trim(bold)//"total k-points: "//&
    trim(reg) , waveFile%kNO
  write(0,"(A,I4)") &
    trim(bold)//"total bands: "//&
    trim(reg) , waveFile%bNO
  write(0,"(A,3I4)") &
    trim(bold)//"nbMAX:"//&
    trim(reg) , (waveFile%nbMAX(j),j=1,3)
  write(0,"(A,3I7)") &
    trim(bold)//"npMAX:"//&
    trim(reg) , (waveFile%npMAX(j),j=1,3)
  write(0,"(A,I7)") &
    trim(bold)//"max plane-waves:"//&
    trim(reg) , waveFile%npCAP
  write(0,"(A,3F11.6)") &
    trim(bold)//"a1 (A):"//&
    trim(reg) , (waveFile%a1(j),j=1,3)
  write(0,"(A,3F11.6)") &
    trim(bold)//"a2 (A):"//&
    trim(reg) , (waveFile%a2(j),j=1,3)
  write(0,"(A,3F11.6)") &
    trim(bold)//"a3 (A):"//&
    trim(reg) , (waveFile%a3(j),j=1,3)
  write(0,"(A,F6.2,A)") &
    trim(bold)//"volume: "//&
    trim(reg) , waveFile%volume , " A^3"
  write(0,"(A,3F11.6)") &
    trim(bold)//"b1 (A):"//&
    trim(reg) , (waveFile%b1(j),j=1,3)
  write(0,"(A,3F11.6)") &
    trim(bold)//"b2 (A):"//&
    trim(reg) , (waveFile%b2(j),j=1,3)
  write(0,"(A,3F11.6)") &
    trim(bold)//"b3 (A):"//&
    trim(reg) , (waveFile%b3(j),j=1,3)
  write(0,"(A,I2)") &
    trim(bold)//"ISPIN: "//&
    trim(reg) , waveFile%spin
  write(0,"(A,I5)") &
    trim(bold)//"precision: "//&
    trim(reg) , waveFile%precision
  write(0,"(A,I6)") &
    trim(bold)//"record length: "//&
    trim(reg) , waveFile%recordLength
  ! reset
  write(0,"(A)") reset
end SUBROUTINE printWAVECAR








! --- WAVECAR pull metadata ------------------------------------------ !
SUBROUTINE pullMetaWAVECAR( waveFile )
  ! workspace
  IMPLICIT none
  TYPE(wavecarFile),INTENT(inout) :: waveFile
  REAL*8,DIMENSION(3) :: b1 , b2 , b3 , TMP
  INTEGER,DIMENSION(3) :: ggA , ggB , ggC
  REAL*8 :: enCut , b1mag , b2mag , b3mag , alpha , beta , kLIMIT
  REAL*8 :: recordIDX , spinIDX , precisionIDX , kIDX , bIDX
  REAL*8 :: xnplane
  INTEGER :: iost , j
  ! set initial values
  waveFile%recordLength = 24
  if (waveFile%unitNo.lt.9) then
    write(0,*) "*** input wavefile unit was below 9 set to 9."
    waveFile%unitNo = 9
  else if (waveFile%unitNo.gt.99) then
    write(0,*) "*** input wavefile unit was above 99 set to 99."
    waveFile%unitNo = 99
  endif
  ! check the file exists
  open( unit=waveFile%unitNo , &
        file=waveFile%filename , &
        access='direct' , &
        recl=waveFile%recordLength , &
        iostat=iost , &
        status='old' )
  if (iost.ne.0) then
    write(6,*) 'error opening file'
    write(0,*) 'Fail opening WAVECAR named: ',trim(waveFile%fileName)
    write(0,*) 'Error code: ' , iost
    write(0,*) 'https://www.ibm.com/support/&
      knowledgecenter/en/SS3KZ4_9.0.0/com.ibm.xlf111.bg.doc/&
      xlflr/iostatvalues.htm'
    close(unit=waveFile%unitNo)
    stop
  endif
  ! pull header variables
  read(unit=waveFile%unitNo,rec=1) recordIDX , spinIDX , precisionIDX
  waveFile%recordLength = nint(recordIDX)
  waveFile%spin = nint(spinIDX)
  waveFile%precision = nint(precisionIDX)
  close(unit=waveFile%unitNo)
  ! only continue if correct precision identified
  if(waveFile%precision.eq.45210) then
    write(0,*) '*** error - WAVECAR_double requires complex*16'
    stop
  endif
  ! re-open file
  open( unit=waveFile%unitNo , &
        file=waveFile%filename , &
        access='direct' , &
        recl=waveFile%recordLength , &
        iostat=iost , &
        status='old' )
  if (iost.ne.0) then
    write(6,*) 'error opening file, on the second go ...'
    write(0,*) 'Failed to open WAVECAR named: ', waveFile%fileName
    write(0,*) 'Error code: ' , iost
    write(0,*) 'https://www.ibm.com/support/&
      knowledgecenter/en/SS3KZ4_9.0.0/com.ibm.xlf111.bg.doc/&
      xlflr/iostatvalues.htm'
    close(unit=waveFile%unitNo)
    stop
  endif
  ! pull vasp variables (kpts,bands,encut,lattice)
  read(unit=waveFile%unitNo,rec=2) kIDX , bIDX , enCut , &
    (waveFile%a1(j),j=1,3) , &
    (waveFile%a2(j),j=1,3) , &
    (waveFile%a3(j),j=1,3)
  waveFile%enCut = nint(enCut)
  waveFile%kNO = nint(kIDX)
  waveFile%bNO = nint(bIDX)
  ! reciprocal lattice
  CALL vcross( waveFile%a2 , waveFile%a3 , b1 )
  CALL vcross( waveFile%a3 , waveFile%a1 , b2 )
  CALL vcross( waveFile%a1 , waveFile%a2 , b3 )
  waveFile%volume = dot_product(waveFile%a1,b1)
  b1 = (2.*pi/waveFile%volume) * b1
  b2 = (2.*pi/waveFile%volume) * b2
  b3 = (2.*pi/waveFile%volume) * b3
  waveFile%b1 = b1
  waveFile%b2 = b2
  waveFile%b3 = b3
  b1mag = dsqrt(dot_product(b1,b1))
  b2mag = dsqrt(dot_product(b2,b2))
  b3mag = dsqrt(dot_product(b3,b3))
  kLIMIT = dsqrt( enCut * c )
  ! A-vec cutoff
  alpha = abs(sin(acos( dot_product(b1,b2) / b1mag / b2mag )))
  CALL vcross(b1,b2 , TMP )
  beta = abs( dot_product(TMP,b3) &
    / dsqrt(dot_product(TMP,TMP)) / b3mag )
  ggA(1) = ( kLIMIT / b1mag / alpha ) + 1
  ggA(2) = ( kLIMIT / b2mag / alpha ) + 1
  ggA(3) = ( kLIMIT / b3mag / beta  ) + 1
  ! B-vec cutoff
  alpha = abs(sin(acos( dot_product(b1,b3) / b1mag / b3mag )))
  CALL vcross(b1,b3 , TMP )
  beta = abs( dot_product(TMP,b2) &
    / dsqrt(dot_product(TMP,TMP)) / b2mag )
  ggB(1) = ( kLIMIT / b1mag / alpha ) + 1
  ggB(2) = ( kLIMIT / b2mag / beta  ) + 1
  ggB(3) = ( kLIMIT / b3mag / alpha ) + 1
  ! C-vec cutoff
  alpha = abs(sin(acos( dot_product(b2,b3) / b2mag / b3mag )))
  CALL vcross(b2,b3 , TMP )
  beta = abs( dot_product(TMP,b1) &
    / dsqrt(dot_product(TMP,TMP)) / b1mag )
  ggC(1) = ( kLIMIT / b1mag / beta  ) + 1
  ggC(2) = ( kLIMIT / b2mag / alpha ) + 1
  ggC(3) = ( kLIMIT / b3mag / alpha ) + 1
  ! formalize cutoffs
  waveFile%nbMAX(1) = max0( ggA(1) , ggB(1) , ggC(1))
  waveFile%nbMAX(2) = max0( ggA(2) , ggB(2) , ggC(2))
  waveFile%nbMAX(3) = max0( ggA(3) , ggB(3) , ggC(3))
  waveFile%npMAX(1) = nint( (4.*pi/3.) * ggA(1)*ggA(2)*ggA(3) )
  waveFile%npMAX(2) = nint( (4.*pi/3.) * ggB(1)*ggB(2)*ggB(3) )
  waveFile%npMAX(3) = nint( (4.*pi/3.) * ggC(1)*ggC(2)*ggC(3) )
  ! set initial planeWave cap
  waveFile%npCAP = min0( &
    waveFile%npMAX(1) , waveFile%npMAX(2) , waveFile%npMAX(3) )
  ! determine magnetic form of WAVECAR
  read( unit=waveFile%unitNo , rec=3 ) xnplane
  if (waveFile%spin.eq.2) then
    waveFile%magState = "collinear"
  else if (nint(xnplane).lt.waveFile%npCAP) then
    waveFile%magState = "nonmagnetic"
  else
    waveFile%magState = "noncollinear"
  endif
  ! done ...
  return
end SUBROUTINE pullMetaWAVECAR







! --- WAVECAR load file from user ------------------------------------ !
SUBROUTINE loadWAVECAR( userArgs , waveFile )
  ! workspace
  IMPLICIT none
  TYPE(userArgumentStruct),INTENT(inout) :: userArgs
  TYPE(wavecarFile),INTENT(inout) :: waveFile
  ! copy file path
  waveFile%fileName = userArgs%fileName_INPUT
  ! load wavecar
  CALL pullMetaWAVECAR( waveFile )
  ! run checks on user selections b,k
  if (userArgs%k.gt.waveFile%kNO) then
    write(0,*) '*** error - selected &
      k=',userArgs%k,' > max k=',waveFile%kNO
    stop
  endif
  if (userArgs%b.gt.waveFile%bNO) then
    write(0,*) '*** error - selected &
      band=',userArgs%b,' > max band=',waveFile%bNO
    stop
  endif
  ! adjust b,k ranges to wavecar
  if (userArgs%kMAX.gt.waveFile%kNO) then
    userArgs%kMAX = waveFile%kNO
  endif
  if (userArgs%bMAX.gt.waveFile%bNO) then
    userArgs%bMAX = waveFile%bNO
  endif
  ! done ...
  return
end SUBROUTINE loadWAVECAR











! --- Acquire raw wavecar data --------------------------------------- !
SUBROUTINE loadWAVEFUNCTION( waveFile , Kpt , Bpt , waveFcn )
  ! --- 0 --- workspace
  IMPLICIT none
  TYPE(wavecarFile),INTENT(in) :: waveFile
  TYPE(waveFunction),INTENT(inout) :: waveFcn
  INTEGER,INTENT(in) :: Kpt , Bpt
  INTEGER :: ncnt , j , ii1,jj1,ii2,jj2,ii3,jj3
  INTEGER :: i , nb1 , nb2 , nb3 , bNO
  REAL*8 :: xnplane
  REAL*8,DIMENSION(3) :: kk
  ! extract variables
  bNO = waveFile%bNO
  nb1 = waveFile%nbMAX(1)
  nb2 = waveFile%nbMAX(2)
  nb3 = waveFile%nbMAX(3)
  ! allocate memory
  if (waveFile%magState.eq.'nonmagnetic') then
    allocate ( waveFcn%uG(waveFile%npCAP,1) )
  else
    allocate ( waveFcn%uG(waveFile%npCAP,2) )
  endif
  allocate ( waveFcn%idx(waveFile%npCAP,3) )
  ! read wavefunction meta-data from wavecar
  read( unit = waveFile%unitNo , rec = 3+(Kpt-1)*(bNO+1) ) &
    xnplane , ( waveFcn%kvec(i) , i=1,3 )
  waveFcn%nplane = nint(xnplane)
  ncnt = 0
  ! LOOP
  do ii3 = 0,2*nb3
  do ii2 = 0,2*nb2
  do ii1 = 0,2*nb1
    ! format planewave index
    jj3 = ii3
    jj2 = ii2
    jj1 = ii1
    if (jj3.gt.nb3) jj3 = jj3 - 2*nb3 - 1
    if (jj2.gt.nb2) jj2 = jj2 - 2*nb2 - 1
    if (jj1.gt.nb1) jj1 = jj1 - 2*nb1 - 1
    ! check energy condition
    do j = 1,3
      kk(j) = &
        (waveFcn%kvec(1)+jj1) * waveFile%b1(j) + &
        (waveFcn%kvec(2)+jj2) * waveFile%b2(j) + &
        (waveFcn%kvec(3)+jj3) * waveFile%b3(j)
    enddo
    ! include planewave index for bounded waves
    if (dot_product(kk,kk)/c.lt.waveFile%encut) then
      ncnt = ncnt+1
      waveFcn%idx(ncnt,1) = jj1
      waveFcn%idx(ncnt,2) = jj2
      waveFcn%idx(ncnt,3) = jj3
    end if
  enddo
  enddo
  enddo
  ! check we counted the appropriate number of planes
  ! NON-COLLINEAR ONLY FOR NOW
  if (2*ncnt.ne.waveFcn%nplane) then
    write(0,*) '*** error - computed 2*ncnt=',2*ncnt, &
        ' != input nplane=',waveFcn%nplane
    write(0,*) 'only noncollinear magnetic wavefunctions for now...'
    stop
  endif
  ! read data
  read( unit=waveFile%unitNo , rec=3+(Kpt-1)*(bNO+1)+Bpt ) &
    ( waveFcn%uG(i,1), i=1,ncnt ) , ( waveFcn%uG(i,2), i=1,ncnt )
  ! save to wavefunction
  waveFcn%nplane = waveFcn%nplane/2
  waveFcn%kpoint = Kpt
  waveFcn%band = Bpt
  waveFcn%src = waveFile
  ! done...
  return
end SUBROUTINE loadWAVEFUNCTION


















! --- help message - user arguments ---------------------------------- !
SUBROUTINE userHelpMessage()
  ! error header
  write(0,*) ''//achar(27)//"[1;3;38;5;088m"//&
    "ERROR reading user commandline arguments - Uh Oh!"
  write(0,*) "Flag syntax (following the executable)..."//&
    achar(27)//"[0;38;5;088m"
  write(0,*) "-f fileName_INPUT   -k kpoint   -K kMIN-kMAX"
  write(0,*) "-b band   -B bMIN-bMAX   -c   -q"
  write(0,*) "-r [fileName_RPOINTS]"
  write(0,*) " "
  write(0,*) ''//achar(27)//"[1;3;38;5;088m"//"Defaults:"//&
              achar(27)//"[0;38;5;088m"
  write(0,*) "fileName_INPUT = WAVECAR"
  write(0,*) "kpoint = 1"
  write(0,*) "kMIN = 1"
  write(0,*) "kMAX = 10000"
  write(0,*) "band = 1"
  write(0,*) "bMIN = 1"
  write(0,*) "bMAX = 10000"
  write(0,*) "fileName_RPOINTS = RPOINTS"
  write(0,*) " "
  write(0,*) ''//achar(27)//"[1;3;38;5;088m"//"NOTE 1: "
  write(0,*) ''//achar(27)//"[0;38;5;088m",&
              "Flags are case-sensitive."
  write(0,*) ''//achar(27)//"[1;3;38;5;088m"//"NOTE 2: "
  write(0,*) ''//achar(27)//"[0;38;5;088m",&
              "Ranges -K and -B can accept '3-' as an example."
  write(0,*) "Make sure the range is in a logical order!"
  write(0,*) ''//achar(27)//"[1;3;38;5;088m"//"NOTE 3: "
  write(0,*) ''//achar(27)//"[0;38;5;088m",&
              "-c Specifies writing 'cartesian' RPOINTS file."
  write(0,*) "If there is an input RPOINTS file, then another one &
              will be written if the -c flag disagrees with the &
              RPOINTS input file format."
  write(0,*) ''//achar(27)//"[1;3;38;5;088m"//"NOTE 4: "
  write(0,*) ''//achar(27)//"[0;38;5;088m",&
              "-q quiets the script, only text output are errors."
  stop
end SUBROUTINE userHelpMessage



! --- file naming (specific to kPoint and band) ---------------------- !
SUBROUTINE makeSuffix( waveFile , k , b , suffix )
  ! --- 0 --- workspace
  IMPLICIT none
  CHARACTER*72,INTENT(out) :: suffix
  INTEGER,INTENT(in) :: k , b
  TYPE(wavecarFile),INTENT(in) :: waveFile
  CHARACTER*72 :: kSTR , kPAD , bSTR , bPAD
  INTEGER :: offset
  ! get values
  write(kSTR,"(I5)") k
  write(bSTR,"(I5)") b
  ! adjust padding
  offset = int(log10(real(waveFile%kNO)))-int(log10(real(k)))
  if (offset.lt.1) then
    kPAD = ""
  else if (offset.lt.2) then
    kPAD = "0"
  else if (offset.lt.3) then
    kPAD = "00"
  else
    kPAD = "000"
  end if
  offset = int(log10(real(waveFile%bNO)))-int(log10(real(b)))
  if (offset.lt.1) then
    bPAD = ""
  else if (offset.lt.2) then
    bPAD = "0"
  else if (offset.lt.3) then
    bPAD = "00"
  else
    bPAD = "000"
  end if
  ! write suffix
  suffix = "_k"//trim(adjustl(kPAD))//trim(adjustl(kSTR))//&
           "_b"//trim(adjustl(bPAD))//trim(adjustl(bSTR))
  ! done...
  return
end SUBROUTINE makeSuffix









end MODULE vtylib__WAVECAR
