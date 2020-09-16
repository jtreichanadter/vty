#! /usr/bin/env python



# ############################## MAIN EXECUTION ############################## #





def parseIncar( *Directory ):

    # specify file-name from Directory
    try:
        Directory = Directory[0]
        fileName = generateFileName(Directory)
    except:
        fileName = generateFileName('.')

    # creates dictionary for storing 'incar' file details
    incar = { key:None for key in keyListIncar() }

    # ##### START LOOP ##### #
    try:
        with open(fileName,'r') as incarFile:
            for line in incarFile:

                # (1) depricate comments
                try: line = line[ :line.index("#") ]
                except ValueError: pass
                line = line.strip()

                # (2) skip over empty (or commented) lines
                if line=="": continue

                # (3) skip lines lines with incorrect '=' count (and report)
                if not line.count("=")==1:
                    warning("Is this INCAR line okay?",line,"(ignoring line)")
                    continue

                # (4a) isolate the 'key'
                key = line[ :line.index("=") ]
                key = key.strip()

                # (4b) isolate the 'value'
                value = line[ line.index("=")+1: ]
                value = value.strip()

                # (5) update INCAR dictionary (or report rogue key)
                if key in keyListIncar():
                    incar[ key ] = value
                else:
                    warning("This INCAR flag is not recognized:",line,"(ignoring line)")

        # post-process raw INCAR dictionary input
        # NOTE: THIS FUNCTION MAY NEED ACTIVE DEVELOPEMENT
        incar = processINCAR(incar)

    except IOError: warning("Evidently there is no INCAR in directory '{0}'".format(Directory))

    # DONE
    return incar










# ############################# HELPER-FUNCTIONS ############################# #


# acquire fileName
def generateFileName(Directory):
    # strip off tail
    Directory = Directory.rstrip('/')
    # base file name
    NAME = "INCAR"
    # if a file, keep as is
    if NAME in Directory:
        fileName = Directory
    # otherwise append file name
    else:
        fileName = Directory + "/" + NAME
    # check if the file exists
    try:
        checkFile = open(fileName,'r')
        checkFile.close()
        return fileName
    except:
        warning("Could not open file '{0}'".format(fileName))
    # try other names
    for j in range(10):
        try:
            checkFile = open(fileName+"-"+str(j),'r')
            checkFile.close()
            return fileName+"-"+str(j)
        except:
            pass
    # could not find the file
    warning("Failed to find {0} file in parse{0}.py!!".format(NAME))
    return NAME


# checks if the given dictionary is generated from an INCAR file
def isDictionaryIncar( dictionary ):
    return \
        not any( not key in keyListIncar() for key in dictionary.keys() ) \
        and \
        not any( not key in dictionary.keys() for key in keyListIncar() )


# prints a full Incar dictionary vertically
def printFullIncar( dictionary ):

    # check that we are looking at an INCAR dictionary
    assert isDictionaryIncar( dictionary )

    # prints everything (including empty keys)
    for item in dictionary.items():
        print("{0}:  {1}".format(item[0],item[1])  )


# prints the nonempty elements of an INCAR dictionary vertically
def printIncar( dictionary ):

    # check that we are looking at an INCAR dictionary
    assert isDictionaryIncar( dictionary )

    # prints everything (including empty keys)
    for item in dictionary.items():
        if not item[1]==None:
            print("{0}: {1}".format(item[0],item[1])  )


# prints the element types of an POSCAR dictionary
def printBriefIncar( dictionary ):
    # check that we are looking at a POSCAR dictionary
    assert isDictionaryIncar( dictionary )
    # prints everything (including empty keys)
    for item in dictionary.items():
        if item[1]==None:
            pass
        elif not isinstance(item[1],list):
            print( "{0}:  {1}".format(item[0],item[1]) )
        else:
            if not isinstance(item[1][0],list):
                print( "{0}:  [{1} list]".format(item[0],len(item[1])) )
            elif not isinstance(item[1][0][0],list):
                print( "{0}:  [{1}x{2} list]".format(item[0],len(item[1]),len(item[1][0])) )
            elif not isinstance(item[1][0][0][0],list):
                print( "{0}:  [{1}x{2}x{3} list]".format(item[0],len(item[1]),len(item[1][0]),len(item[1][0][0])) )
            elif not isinstance(item[1][0][0][0][0],list):
                print( "{0}:  [{1}x{2}x{3}x{4} list]".format(item[0],len(item[1]),len(item[1][0]),len(item[1][0][0]),len(item[1][0][0][0])) )
            else:
                print( "{0}:  WHAT IS THIS??".format(item[0]) )


# gets type of INCAR run (and prints result)
def getTypeIncar( dictionary ):
    # This is currently primative as it checks just for a couple of flags,
    # not entriely exhaustive.
    # Note: dictionary must be from the parseIncar() function.
    #
    # returns
    # 0 - error!
    # 1 - scf
    # 2 - relax
    # 4 - bandstructure

    # check for error (must have matching keys with the INCAR flag list)
    if not isDictionaryIncar( dictionary ):
        print("\033[1;33m--- WARNING : getTypeIncar() ---\033[0;33m")
        print("This doesn't appear to be a valid INCAR dictionary")
        print("-- May not have been accessed by processINCAR()?")
        printDictionary( dictionary )
        print("\033[0m")
        return 0

    # check IBRION key for relaxation
    if not dictionary["IBRION"]==None and dictionary["IBRION"]!=-1:
        print("\033[1;37mThis INCAR is for a relaxation!\033[0m")
        return 2

    # check ICHARGE key for bandstructures
    if not dictionary["ICHARG"]==None and dictionary["ICHARG"]>9:
        print("\033[1;37mThis INCAR is for a bandstructure calculation!\033[0m")
        return 4

    # by elemination it's an scf calculation
    print("\033[1;37mThis INCAR is for an scf calculation!\033[0m")
    return 1



# prints INCAR parser warnings
def warning(*argv):
    print("\033[1;33m--- WARNING : parseIncar ---\033[0;33m")
    for arg in argv:
        print(arg)
    print("\033[0m")
















# ############################## INCAR-HELPERS ############################### #



# returns an alphabetical list of incar keys
def keyListIncar():

    # ### INITIALIZE DICTIONARY
    # These keys are based on the online VASP wiki documentation.
    # --> URL: https://www.vasp.at/wiki/index.php/Category:INCAR
    # --> "This page was last edited on 22 February 2019"
    # --> Each row contains 5 tags; for 345 total entries.
    # --> Last updated: December 3rd, 2019

    return [ "ADDGRID" , "AEXX" , "AGGAC" , "AGGAX" , "ALDAC" , \
    "ALGO" , "AMIN" , "AMIX" , "AMIX MAG" , "ANDERSEN PROB" , \
    "ANTIRES" , "APACO" , "BMIX" , "BMIX MAG" , "CH LSPEC" , \
    "CH NEDOS" , "CH SIGMA" , "CLL" , "CLN" , "CLNT" , \
    "CLZ" , "CMBJ" , "CMBJA" , "CMBJB" , "CSHIFT" , \
    "DEPER" , "DIMER DIST" , "DIPOL" , "DQ" , \
    "EBREAK", "EDIFF", "EDIFFG", "EFIELD", "EFIELD PEAD" , \
    "EINT" , "Electric Field Gradient" , "EMAX" , "EMIN" , "ENAUG" , \
    "ENCUT" , "ENCUTFOCK" , "ENCUTGW" , "ENCUTGWSOFT" , "ENINI" , \
    "EPSILON" , "EVENONLY" , "EVENONLYGW" , "FERDO" , "FERWE" , \
    "FINDIFF" , "GGA" , "GGA COMPAT" , "HFLMAX" , "HFRCUT" , \
    "HFSCREEN" , "HILLS BIN" , "HILLS H" , "HILLS W" , "HITOLER" \
    "I CONSTRAINED M" , "IALGO" ,"IBAND" , "IBRION" , "ICHARG" , \
    "ICHIBARE" , "ICORELEVEL" , "IDIPOL" , "IEPSILON" , "IGPAR" , \
    "IMAGES" , "IMIX" , "INCREM" , "INIMIX" , "INIWAV" , \
    "IPEAD" , "ISIF" , "ISMEAR" , "ISPIN" , "ISTART" , \
    "ISYM" , "IVDW" , "IWAVPR" , "KBLOCK" , "KGAMMA" , \
    "KPAR" , "KPOINT BSE" , "KPUSE" , "KSPACING" , "LADDER" , \
    "LAECHG" , "LAMBDA" , "LANGEVIN GAMMA" , "LANGEVIN GAMMA L" , "LASPH" , \
    "LASYNC" , "LATTICE CONSTRAINTS" , "LBERRY" , "LBLUEOUT" , "LBONE" , \
    "LCALCEPS" , "LCALCPOL" , "LCHARG" , "LCHIMAG" , "LCORR" , \
    "LDAU" , "LDAUJ" , "LDAUL" , "LDAUPRINT" , "LDAUTYPE" , \
    "LDAUU" , "LDIAG" , "LDIPOL" , "LEFG" , "LELF" , \
    "LEPSILON" , "LFOCKAEDFT" , "LHARTREE" , "LHFCALC" , "LHYPERFINE" , \
    "LKPROJ" , "LLRAUG" , "LMAXFOCK" , "LMAXFOCKAE" , "LMAXMIX" , \
    "LMAXPAW" , "LMAXTAU" , "LMIXTAU" , "LMONO" , "LNABLA" , \
    "LNMR SYM RED" , "LNONCOLLINEAR" , "LOCPROJ" , "LOPTICS" , "LORBIT" , \
    "LORBMOM" , "LPARD" , "LPEAD" , "LPLANE" , "LREAL" , \
    "LRPA" , "LSCAAWARE" , "LSCALAPACK" , "LSCALU" , "LSCSGRAD" , \
    "LSELFENERGY" , "LSEPB" , "LSEPK" , "LSORBIT" , "LSPECTRAL" , \
    "LSPECTRALGW" , "LSPIRAL" , "LSUBROT" , "LTHOMAS" , "LUSE VDW" , \
    "LVDW EWALD" , "LVDW ONECELL" , "LVDWEXPANSION" , "LVHAR" , "LVTOT" , \
    "LWANNIER90" , "LWANNIER90 RUN" , "LWAVE" , "LWRITE MMN AMN" , "LWRITE UNK" , \
    "LWRITE WANPROJ" , "LZEROZ" , "M CONSTR" , "MAGMOM" , "MAXMEM" , "MAXMIX" , \
    "MDALGO" , "METAGGA" , "MINROT" , "MIXPRE" , "ML FF AFILT2 MB" , \
    "ML FF CDOUB" , "ML FF CSF" , "ML FF CSIG" , "ML FF CSLOPE" , "ML FF CTIFOR" , \
    "ML FF EATOM" , "ML FF IAFILT2 MB" , "ML FF IBROAD1 MB" , "ML FF IBROAD2 MB" , "ML FF ICOUPLE MB" , \
    "ML FF ICUT1 MB" , "ML FF ICUT2 MB" , "ML FF IERR" , "ML FF IREG MB" , "ML FF ISAMPLE" , \
    "ML FF ISCALE TOTEN MB" , "ML FF ISOAP1 MB" , "ML FF ISOAP2 MB" , "ML FF ISTART" , "ML FF IWEIGHT" , \
    "ML FF LAFILT2 MB" , "ML FF LBASIS DISCARD" , "ML FF LCONF DISCARD" , "ML FF LCOUPLE MB" , "ML FF LCRITERIA" , \
    "ML FF LEATOM MB" , "ML FF LHEAT MB" , "ML FF LMAX2 MB" , "ML FF LMLFF" , "ML FF LMLMB" , \
    "ML FF LNORM1 MB" , "ML FF LNORM2 MB" , "ML FF MB MB" , "ML FF MCONF" , "ML FF MCONF NEW" , \
    "ML FF MHIS" , "ML FF MRB1 MB" , "ML FF MRB2 MB" , "ML FF MSPL1 MB" , "ML FF MSPL2 MB" , \
    "ML FF NATOM COUPLED MB" , "ML FF NDIM SCALAPACK" , "ML FF NHYP1 MB" , "ML FF NHYP2 MB" , "ML FF NMDINT" , \
    "ML FF NR1 MB" , "ML FF NR2 MB" , "ML FF NWRITE" , "ML FF RCOUPLE MB" , "ML FF RCUT1 MB" , \
    "ML FF RCUT2 MB" , "ML FF SIGV0 MB" , "ML FF SIGW0 MB" , "ML FF SION1 MB" , "ML FF SION2 MB" , \
    "ML FF W1 MB" , "ML FF W2 MB" , "ML FF WTIFOR" , "ML FF WTOTEN" , "ML FF WTSIF" , \
    "NBANDS" , "NBANDSGW" , "NBANDSO" , "NBANDSV" , "NBLK" , \
    "NBLOCK" , "NBMOD" , "NBSEEIG" , "NCORE" , "NCRPA BANDS" , \
    "NDAV" , "NEDOS" , "NELECT" , "NELM" , "NELMDL" , \
    "NELMIN" , "NFREE" , "NGX" , "NGXF" , "NGY" , \
    "NGYF" , "NGYROMAG" , "NGZ" , "NGZF" , "NKRED" , \
    "NKREDX" , "NKREDY" , "NKREDZ" , "NLSPLINE" , "NMAXFOCKAE" , \
    "LMAXFOCKAE" , "NOMEGA" , "NOMEGAPAR" , "NOMEGAR" , "NPACO" , \
    "NPAR" , "NPPSTR" , "NSIM" , "NSUBSYS" , "NSW" , \
    "NTARGET STATES" , "NTAUPAR" , "NUPDOWN" , "NWRITE" , "ODDONLY" , \
    "ODDONLYGW" , "OFIELD A" , "OFIELD KAPPA" , "OFIELD Q6 FAR" , "OFIELD Q6 NEAR" , \
    "OMEGAMAX" , "OMEGAMIN" , "OMEGATL" , "PARAM1" , "PARAM2" , \
    "PFLAT" , "PHON LBOSE" , "PHON LMC" , "PHON NSTRUCT" , "PHON NTLIST" \
    "PHON TLIST" , "PLEVEL" , "PMASS" , "POMASS" , "POTIM" , \
    "PREC" , "PRECFOCK" , "PROUTINE" , "PSTRESS" , "PSUBSYS" , \
    "PTHRESHOLD" , "QMAXFOCKAE" , "QSPIRAL" , "QUAD EFG" , "RANDOM SEED" , \
    "ROPT" , "RWIGS" , "SAXIS" , "SCSRAD" , "SHAKEMAXITER" , \
    "SHAKETOL" , "SIGMA" , "SMASS" , "SMEARINGS" , "SPRING" , \
    "STEP MAX" , "STEP SIZE" , "SYMPREC" , "SYSTEM" , "TEBEG" , \
    "TEEND" , "TIME" , "TSUBSYS" , "VALUE MAX" , "VALUE MIN" , \
    "VCUTOFF" , "VDW A1" , "VDW A2" , "VDW C6" , "VDW CNRADIUS" , \
    "VDW D" , "VDW R0" , "VDW RADIUS" , "VDW S6" , "VDW S8" , \
    "VDW SR" , "VOSKOWN" , "WC" , "WEIMIN" , "ZVAL" ]






# converts raw strings in INCAR dictionary (after parsing) to their proper values
# REQUIRES REGULAR DEVELOPMENT FOR NEW KEYS
def processINCAR( dictionary ):
    # ### PROCESS RAW INCAR DICTIONARY
    # Documentation and allowed values stem from the online VASP wiki.
    # --> URL: https://www.vasp.at/wiki/index.php/Category:INCAR
    # --> "This page was last edited on 22 February 2019"
    # --> Last updated: December 3rd, 2019

    # loop through nonempty dictionary values
    for key,value in dictionary.items():
        if not value == None:

            #	ADDGRID
            #	AEXX
            #	AGGAC
            #	AGGAX
            #	ALDAC

            # ALGO - DFT algorithm
            if key=="ALGO":
                options = ["Normal","VeryFast","Fast","Conjugate","All","Damped","Subrot","Eigenval","Exact","None","Nothing","CHI","G0W0","GW0","GW","scGW0","scGW","G0W0R","GW0R","GWR","scGW0R","scGWR","ACFDT","RPA","ACFDTR","RPAR","BSE","TDHF"]
                if not any( value==option for option in options ):
                    warning("Failed to assign '{0}' with value '{1}'.".format(key,value),"Available options: {0}".format(options))
                    dictionary[key] = None
                # Normal: selects IALGO=38 (blocked Davidson iteration scheme).
                # VeryFast: selects IALGO=48 (RMM-DIIS).
                # Fast: selects a faily robust mixture of the Davidson and RMM-DIIS algorithms. In this case, Davidson (IALGO=38) is used for the initial phase, and then VASP switches to RMM-DIIS (IALGO=48). Subsequencly, for each ionic update, one IALGO=38 sweep is performed for each ionic step (except the first one).
                # Conjugate or All: selects an "all band simultaneous update of orbitals" (IALGO=58, in both cases the same conjugate gradient algorithm is used).
                # Damped: selects a damped velocity friction algorithm (IALGO=53).
                # Exact: performs an exact diagonalization (IALGO=90).
                # Subrot: selects subspace rotation or diagonalization in the sub-space spanned by the orbitals (IALGO=4).
                # Eigenval: allows to recalculate one electron energies, density of state and perform selected postprocessing using the current orbitals (IALGO=3) e.g. read from the WAVECAR file.
                # None or Nothing: allows to recalculate the density of states or perform selected postprocessing, using the current orbitals and one electron energies (IALGO=2) e.g. read from the WAVECAR file.
                #   (( Following tags are available as of VASP.5.X. ))
                # CHI: calculates the response functions only.
                # TDHF: selects TDHF calculations using the VASP internal Cassida code see BSE calculations, (available as of VASP.5.2.12)
                # BSE selects BSE calculations using the VASP internal Cassida code (see BSE calculations)
                # Timeev: performs a delta-pulse in time and then performs timepropagation
                # ACFDT: selects RPA total energy calculations see ACFDT/RPA calculations
                # RPA: synonymous to ACFDT see ACFDT/RPA calculations (available as of VASP.5.3.1)
                #   (( available as of VASP.6 ))
                # RPAR: selects low scaling RPA total energy calculations (for details see ACFDT/RPA calculations)
                # ACFTDR: synonym for RPAR (for details see ACFDT/RPA calculations)
                # QPGW0R: selects low scaling version of QPGW0. Quasi-particle GW calculations are performed, where off-diagonal components of the self-energy are included. A full update of the QP energies AND one electron orbitals is performed in the calculation of G only (for details see QPGW0R calculations)
                # QPGWR: selects low scaling version of QPGW. Quasi-particle GW calculations are performed, where off-diagonal components of the self-energy are included. A full update of the QP energies AND one electron orbitals is performed in the calculation of G and W (for details see QPGWR calculations)
                # scGW0R: selects self-consistent GW0 calculations, where only the Green's function G is updated from the corresponding Dyson. The screened potential W remains unchanged after the first iteration. NELM iteration cycles are performed (see self-consistent GW calculations).
                # scGWR: selects self-consistent GW calculations, where both, G and W are updated from the corresponding Dyson equation. NELM iteration cycles are performed. (for details see self-consistent GW calculations).
                continue

            #	AMIN
            #	AMIX
            #	AMIX MAG
            #	ANDERSEN PROB
            #	ANTIRES
            #	APACO
            #	BMIX
            #	BMIX MAG
            #	CH LSPEC
            #	CH NEDOS
            #	CH SIGMA
            #	CLL
            #	CLN
            #	CLNT
            #	CLZ
            #	CMBJ
            #	CMBJA
            #	CMBJB
            #	CSHIFT
            #	DEPER
            #	DIMER DIST
            #	DIPOL
            #	DQ
            #	EBREAK

            # EDIFF - elecronic convergence criteria (eV)
            if key=="EDIFF":
                dictionary[key] = float( value )
                if dictionary[key]<1e-6 or dictionary[key]>1e-2 and not dictionary[key]==0:
                    warning("Did you mean to use this electronic convergence?","EDIFF = {0}".format(value))
                # For EDIFF=0, NELM electronic SC-steps will always be performed
                continue

            # EDIFFG - ionic convergence critera (eV/A)
            if key=="EDIFFG":
                dictionary[key] = float( value )
                if dictionary[key]<1e-6 or dictionary[key]>1e-2 and not dictionary[key]==0:
                    warning("Did you mean to use this ionic convergence?","EDIFFG = {0}".format(value))
                # For EDIFG=0, NSW steps will always be performed
                continue

            #	EFIELD
            #	EFIELD PEAD
            #	EINT
            #	Electric Field Gradient
            #	EMAX
            #	EMIN
            #	ENAUG

            # ENCUT - plane wave basis frequency limit
            if key=="ENCUT":
                dictionary[key] = int(float( value ))
                # Gcut = sqrt( ENCUT * 2m/bb^2 )
                continue

            #	ENCUTFOCK
            #	ENCUTGW
            #	ENCUTGWSOFT
            #	ENINI
            #	EPSILON
            #	EVENONLY
            #	EVENONLYGW
            #	FERDO
            #	FERWE
            #	FINDIFF
            #	GGA
            #	GGA COMPAT
            #	HFLMAX
            #	HFRCUT
            #	HFSCREEN
            #	HILLS BIN
            #	HILLS H
            #	HILLS W
            #	HITOLER
            #	I CONSTRAINED M
            #	IALGO
            #	IBAND

            # IBRION - determines how the ions are updated and moved.
            if key=="IBRION":
                dictionary[key] = int(float( value ))
                options = [-1,0,1,2,3,5,6,7,8,44]
                if not any( dictionary[key]==option for option in options ):
                    warning("Failed to assign '{0}' with value '{1}'.".format(key,value),"Available options: {0}".format(options))
                    dictionary[key] = None
                # IBRION=-1: no update.
                # IBRION=0: molecular dynamics.
                # IBRION=1: ionic relaxation (RMM-DIIS).
                # IBRION=2: ionic relaxation (conjugate gradient algorithm).
                # IBRION=3: ionic relaxation (damped molecular dynamics).
                # IBRION=5 and 6: second derivatives, Hessian matrix and phonon frequencies (finite differences).
                # IBRION=7 and 8: second derivatives, Hessian matrix and phonon frequencies (perturbation theory).
                # IBRION=44: the Improved Dimer Method.
                continue

            # ICHARG - determines how VASP constructs initial charge density.
            if key=="ICHARG":
                dictionary[key] = int(float( value ))
                options = [0,1,2,4,10,11,12,14]
                if not any( dictionary[key]==option for option in options ):
                    warning("Failed to assign '{0}' with value '{1}'.".format(key,value),"Available options: {0}".format(options))
                    dictionary[key] = None
                # 0             Calculate charge density from initial wave functions.
                #                 If ISTART is internally reset due to an invalid WAVECAR-file,
                #                 ICHARG will be set to ICHARG=2.
                # 1             Read the charge density from file CHGCAR,
                #                 extrapolate from the old positions (on CHGCAR) to
                #                 the new positions using a linear combination of
                #                 atomic charge densities.
                # 2             Take superposition of atomic charge densities.
                # 4             Read potential from file POT (VASP 5.1)
                # 10,11,12,14   Adding 10 holds the charge density constant during
                #                 the whole electronic minimization.
                continue

            #	ICHIBARE
            #	ICORELEVEL
            #	IDIPOL
            #	IEPSILON
            #	IGPAR
            #	IMAGES
            #	IMIX
            #	INCREM
            #	INIMIX
            #	INIWAV
            #	IPEAD

            # ISIF - sets constraints on ionic relaxation steps
            if key=="ISIF":
                dictionary[key] = int(float( value ))
                options = [0,1,2,3,4,5,6,7]
                if not any( dictionary[key]==option for option in options ):
                    warning("Failed to assign '{0}' with value '{1}'.".format(key,value),"Available options: {0}".format(options))
                    dictionary[key] = None
                #       calculate                 degrees-of-freedom
                #   forces	Stress tensor	positions	cell shape	cell volume
                # 0	yes     no              yes         no          no
                # 1	yes     trace only      yes         no          no
                # 2	yes     yes             yes         no          no
                # 3	yes     yes             yes         yes         yes
                # 4	yes     yes             yes         yes         no
                # 5	yes     yes             no          yes         no
                # 6	yes     yes             no          yes         yes
                # 7	yes     yes             no          no          yes
                continue

            # ISMEAR - determines the partial occupancies 'fnk' for each orbital.
            if key=="ISMEAR":
                dictionary[key] = int(float( value ))
                options = [-5,-4,-3,-2,-1]
                if not any( dictionary[key]==option for option in options ) and dictionary[key]<0:
                    warning("Failed to assign '{0}' with value '{1}'.".format(key,value),"Available options: {0}".format(options))
                    dictionary[key] = None
                # ISMEAR=N (N>0): method of Methfessel-Paxton order N (For the Methfessel-Paxton scheme the partial occupancies can be negative)
                # ISMEAR=0: Gaussian smearing (USE FOR BANDSTRUCTURE)
                # ISMEAR=−1: Fermi smearing.
                # ISMEAR=−2: partial occupancies are read in from the WAVECAR or INCAR file, and kept fixed throughout run. (occupancies also written to OUTCAR; in this case they are multiplied by 2, i.e. between 0 and 2)
                # ISMEAR=−3: perform a loop over smearing-parameters supplied in the INCAR file. MUST WRITE 'SMEARINGS' TAG.
                # ISMEAR=−4: tetrahedron method (use a Γ-centered k-mesh).
                # ISMEAR=−5: tetrahedron method with Blöchl corrections (use a Γ-centered k-mesh).
                # (note that SIGMA determines the width of the smearing, eV)
                continue

            # ISPIN - specifies spin polarization
            if key=="ISPIN":
                dictionary[key] = int(float( value ))
                if not ( dictionary[key]==1 or dictionary[key]==2 ):
                    warning("Failed to assign '{0}' with value '{1}'.".format(key,value),"Available options are: 1,2")
                    dictionary[key] = None
                # ISPIN=1: non spin polarized calculations are performed
                # ISPIN=2: spin polarized calculations (collinear) are performed
                # (can study collinear magnetism by combining ISPIN with MAGMOM )
                continue

            #	ISTART

            # ISYM - determines how VASP handles symmetry
            if key=="ISYM":
                dictionary[key] = int(float( value ))
                options = [0,1,2,3,4,5,6,7]
                if not any( dictionary[key]==option for option in options ):
                    warning("Failed to assign '{0}' with value '{1}'.".format(key,value),"Available options: {0}".format(options))
                    dictionary[key] =  None
                # -1 the use of symmetry is switched off completely.
                #  0, VASP does not use symmetry, but it will assume that Ψk=Ψ*-k and reduces the sampling of the Brillouin zone.
                #  1, symmetrization is on.
                #  2, a more efficient, memory conserving symmetrisation of the charge density is used.
                #  3, VASP does not directly symmetrize the charge density. Instead the charge density is constructed by applying the relevant symmetry operations to the orbitals at the k-points in the irreducible part of the Brillouin zone. This method of symmetrization is used when LHFCALC=.TRUE
                continue

            # IVDW - whether vdW corrections are calculated or not (and how).
            if key=="IVDW":
                dictionary[key] = int(float( value ))
                options = [0,1,10,11,12,2,20,21,202,4]
                if not any( dictionary[key]==option for option in options ):
                    warning("Failed to assign '{0}' with value '{1}'.".format(key,value),"Available options: {0}".format(options))
                    dictionary[key] = None
                # 0     no correction
                # 1|10  DFT-D2 method of Grimme (available as of VASP.5.2.11)
                # 11    zero damping DFT-D3 method of Grimme (available as of VASP.5.3.4)
                # 12    DFT-D3 method with Becke-Jonson damping (available as of VASP.5.3.4)
                # 2|20  Tkatchenko-Scheffler method (available as of VASP.5.3.3)
                # 21    Tkatchenko-Scheffler method with iterative Hirshfeld partitioning (available as of VASP.5.3.5)
                # 202   Many-body dispersion energy method (MBD@rSC) (available as of VASP.5.4.1)
                # 4     dDsC dispersion correction method (available as of VASP.5.4.1)
                continue

            #	IWAVPR
            #	KBLOCK
            #	KGAMMA
            #	KPAR
            #	KPOINT BSE
            #	KPUSE
            #	KSPACING
            #	LADDER
            #	LAECHG
            #	LAMBDA
            #	LANGEVIN GAMMA
            #	LANGEVIN GAMMA L
            #	LASPH
            #	LASYNC
            #	LATTICE CONSTRAINTS
            #	LBERRY
            #	LBLUEOUT
            #	LBONE
            #	LCALCEPS
            #	LCALCPOL
            #	LCHARG
            #	LCHIMAG
            #	LCORR
            #	LDAU
            #	LDAUJ
            #	LDAUL
            #	LDAUPRINT
            #	LDAUTYPE
            #	LDAUU
            #	LDIAG
            #	LDIPOL
            #	LEFG
            #	LELF
            #	LEPSILON
            #	LFOCKAEDFT
            #	LHARTREE
            #	LHFCALC
            #	LHYPERFINE
            #	LKPROJ
            #	LLRAUG
            #	LMAXFOCK
            #	LMAXFOCKAE
            #	LMAXMIX
            #	LMAXPAW
            #	LMAXTAU
            #	LMIXTAU
            #	LMONO
            #	LNABLA
            #	LNMR SYM RED

            # LNONCOLLINEAR - whether fully non-collinear magnetic calculations are performed
            if key=="LNONCOLLINEAR":
                if   value==".TRUE.": dictionary[key] = True
                elif value==".FALSE": dictionary[key] = False
                else:
                    warning("Failed to assign '{0}' with value '{1}'.".format(key,value),"Available options: .TRUE. .FALSE.")
                    dictionary[key] = None
                # VASP is capable of reading WAVECAR and CHGCAR files from previous non-magnetic or collinear calculations
                # For a non-collinear setup, three values must be supplied for each ion in the MAGMOM line
                continue

            #	LOCPROJ
            #	LOPTICS

            # LORBIT - together with an appropriate RWIGS, determines whether the PROCAR or PROOUT files are written
            if key=="LORBIT":
                dictionary[key] = int(float( value ))
                options = [0,1,2,5,10,11,12]
                if not any( dictionary[key]==option for option in options ):
                    warning("Failed to assign '{0}' with value '{1}'.".format(key,value),"Available options: {0}".format(options))
                    dictionary[key] = None
                # LORBIT	RWIGS tag	files written
                # 0         required	DOSCAR & PROCAR
                # 1         required	DOSCAR & lm-decomposed PROCAR
                # 2         required	DOSCAR & lm-decomposed PROCAR + phase factors
                # 5         required	DOSCAR & PROOUT
                # 10        ignored     DOSCAR & PROCAR
                # 11        ignored     DOSCAR & lm-decomposed PROCAR
                # 12        ignored     DOSCAR & lm-decomposed PROCAR + phase factors
                continue

            #	LORBMOM
            #	LPARD
            #	LPEAD

            # LPLANE - switches on the plane-wise data distribution in real space
            if key=="LPLANE":
                if   value==".TRUE.": dictionary[key] = True
                elif value==".FALSE": dictionary[key] = False
                else:
                    warning("Failed to assign '{0}' with value '{1}'.".format(key,value),"Available options: .TRUE. .FALSE.")
                    dictionary[key] = None
                    # general algorithmic flag
                    # .TRUE. reduces the communication band width during the FFT's,
                    #        but worsens the load balancing on massively parallel machines.
                    #        Only use if NGZ is >3×(number of nodes)/NPAR, and optimal load balancing is achieved if NGZ=n×NPAR.
                continue

            #	LREAL
            #	LRPA
            #	LSCAAWARE
            #	LSCALAPACK
            #	LSCALU
            #	LSCSGRAD
            #	LSELFENERGY
            #	LSEPB
            #	LSEPK

            # LSORBIT - specifies whether spin-orbit coupling is taken into account
            if key=="LSORBIT":
                if   value==".TRUE.": dictionary[key] = True
                elif value==".FALSE": dictionary[key] = False
                else:
                    warning("Failed to assign '{0}' with value '{1}'.".format(key,value),"Available options: .TRUE. .FALSE.")
                    dictionary[key] = None
                # Only works for PAW potentials; not supported by ultrasoft pseudopotentials.
                # If spin-orbit coupling is not included, energy is moment direction independent.
                continue

            #	LSPECTRAL
            #	LSPECTRALGW
            #	LSPIRAL
            #	LSUBROT
            #	LTHOMAS
            #	LUSE VDW
            #	LVDW EWALD
            #	LVDW ONECELL
            #	LVDWEXPANSION
            #	LVHAR
            #	LVTOT
            #	LWANNIER90
            #	LWANNIER90 RUN
            #	LWAVE
            #	LWRITE MMN AMN
            #	LWRITE UNK
            #	LWRITE WANPROJ
            #	LZEROZ
            #	M CONSTR

            # MAGMOM - Specifies the initial magnetic moment for each atom.
            if key=="MAGMOM":
                # initialize MAGMOM array
                dictionary[key] = []
                # sweep through listed entries
                entries = value.split()
                for entry in entries:
                    print
                    try:
                        numMags = int(float( entry[ :entry.index("*") ] ))
                        magValue = float( entry[ entry.index("*")+1: ] )
                        for j in range(numMags):
                            dictionary[key].append(magValue)
                    except ValueError: dictionary[key].append(float( entry ))
                # Works ONLY if ICHARG=2, or if ICHARG=1 and the CHGCAR file contains no magnetisation density.
                # Start from larger local magnetic moments.
                #    A safe default is usually the experimental magnetic moment multiplied by 1.2 or 1.5.
                continue

            #	MAXMEM
            #	MAXMIX
            #	MDALGO
            #	METAGGA
            #	MINROT
            #	MIXPRE
            #	ML FF AFILT2 MB
            #	ML FF CDOUB
            #	ML FF CSF
            #	ML FF CSIG
            #	ML FF CSLOPE
            #	ML FF CTIFOR
            #	ML FF EATOM
            #	ML FF IAFILT2 MB
            #	ML FF IBROAD1 MB
            #	ML FF IBROAD2 MB
            #	ML FF ICOUPLE MB
            #	ML FF ICUT1 MB
            #	ML FF ICUT2 MB
            #	ML FF IERR
            #	ML FF IREG MB
            #	ML FF ISAMPLE
            #	ML FF ISCALE TOTEN MB
            #	ML FF ISOAP1 MB
            #	ML FF ISOAP2 MB
            #	ML FF ISTART
            #	ML FF IWEIGHT
            #	ML FF LAFILT2 MB
            #	ML FF LBASIS DISCARD
            #	ML FF LCONF DISCARD
            #	ML FF LCOUPLE MB
            #	ML FF LCRITERIA
            #	ML FF LEATOM MB
            #	ML FF LHEAT MB
            #	ML FF LMAX2 MB
            #	ML FF LMLFF
            #	ML FF LMLMB
            #	ML FF LNORM1 MB
            #	ML FF LNORM2 MB
            #	ML FF MB MB
            #	ML FF MCONF
            #	ML FF MCONF NEW
            #	ML FF MHIS
            #	ML FF MRB1 MB
            #	ML FF MRB2 MB
            #	ML FF MSPL1 MB
            #	ML FF MSPL2 MB
            #	ML FF NATOM COUPLED MB
            #	ML FF NDIM SCALAPACK
            #	ML FF NHYP1 MB
            #	ML FF NHYP2 MB
            #	ML FF NMDINT
            #	ML FF NR1 MB
            #	ML FF NR2 MB
            #	ML FF NWRITE
            #	ML FF RCOUPLE MB
            #	ML FF RCUT1 MB
            #	ML FF RCUT2 MB
            #	ML FF SIGV0 MB
            #	ML FF SIGW0 MB
            #	ML FF SION1 MB
            #	ML FF SION2 MB
            #	ML FF W1 MB
            #	ML FF W2 MB
            #	ML FF WTIFOR
            #	ML FF WTOTEN
            #	ML FF WTSIF
            #	NBANDS
            #	NBANDSGW
            #	NBANDSO
            #	NBANDSV
            #	NBLK
            #	NBLOCK
            #	NBMOD
            #	NBSEEIG

            # NCORE - determines the number of compute cores that work on an individual orbital (available as of VASP.5.2.13)
            if key=="NCORE":
                dictionary[key] = int(float( value ))
                if dictionary[key]<1:
                    warning("Failed to assign '{0}' with value '{1}'.".format(key,value),"Available options: .TRUE. .FALSE.")
                    dictionary[key] = None
                # NCORE = number-of-cores per node OR sqrt(number-of-cores)
                continue

            #	NCRPA BANDS
            #	NDAV
            #	NEDOS
            #	NELECT

            # NELM - sets the maximum number of electronic SCF steps
            if key=="NELM":
                dictionary[key] = int(float( value ))
                if dictionary[key]<1:
                    warning("Failed to assign '{0}' with value '{1}'.".format(key,value),"Available options: .TRUE. .FALSE.")
                    dictionary[key] = None
                # if you're not converging after 40, reconsider IALGO or ALGO, LSUBROT, and the mixing-parameters.
                continue

            #	NELMDL

            #	NELMIN - sets the minimum number of electronic SCF steps
            if key=="NELMIN":
                dictionary[key] = int(float( value ))
                if dictionary[key]<1:
                    warning("Failed to assign '{0}' with value '{1}'.".format(key,value),"Available options: .TRUE. .FALSE.")
                    dictionary[key] = None
                continue

            #	NFREE
            #	NGX
            #	NGXF
            #	NGY
            #	NGYF
            #	NGYROMAG
            #	NGZ
            #	NGZF
            #	NKRED
            #	NKREDX
            #	NKREDY
            #	NKREDZ
            #	NLSPLINE
            #	NMAXFOCKAE
            #	LMAXFOCKAE
            #	NOMEGA
            #	NOMEGAPAR
            #	NOMEGAR
            #	NPACO
            #	NPAR
            #	NPPSTR

            # NSIM - sets number of bands optimized simultaneously by RMM-DIIS algorithm
            if key=="NSIM":
                dictionary[key] = int(float( value ))
                if dictionary[key]<1:
                    warning("Failed to assign '{0}' with value '{1}'.".format(key,value),"Available options: .TRUE. .FALSE.")
                    dictionary[key] = None
                continue

            #	NSUBSYS

            # NSW - sets the maximum number of ionic steps
            if key=="NSW":
                dictionary[key] = int(float( value ))
                # most relevant when IBRION!=0
                continue

            #	NTARGET STATES
            #	NTAUPAR
            #	NUPDOWN

            # NWRITE - determines how much will be written to the file OUTCAR ('verbosity flag').
            if key=="NWRITE":
                dictionary[key] = int(float( value ))
                options = [0,1,2,3,4]
                if not any( dictionary[key]==option for option in options ):
                    warning("Failed to assign '{0}' with value '{1}'.".format(key,dictionary[key]),"Available options: {0}".format(options))
                    dictionary[key] = None
                # Default NWRITE=2
                # For long MD-runs use NWRITE=0 or NWRITE=1.
                # For short runs use NWRITE=2.
                # NWRITE=3 might give information if something goes wrong.
                # NWRITE=4 is for debugging only
                continue

            #	ODDONLY
            #	ODDONLYGW
            #	OFIELD A
            #	OFIELD KAPPA
            #	OFIELD Q6 FAR
            #	OFIELD Q6 NEAR
            #	OMEGAMAX
            #	OMEGAMIN
            #	OMEGATL
            #	PARAM1
            #	PARAM2
            #	PFLAT
            #	PHON LBOSE
            #	PHON LMC
            #	PHON NSTRUCT
            #	PHON NTLIST
            #	PHON TLIST
            #	PLEVEL
            #	PMASS
            #	POMASS
            #	POTIM

            # PREC - specifies the "precision"-mode.
            if key=="PREC":
                options = ["Low","Medium","High","Normal","Single","Accurate"]
                if not any( value==option for option in options ):
                    warning("Failed to assign '{0}' with value '{1}'.".format(key,value),"Available options",options)
                    dictionary[key] = None
                # The PREC parameter sets the defaults for four sets of parameters (ENCUT; NGX, NGY, NGZ; NGXF, NGYF, NGZF, and ROPT)
                # If very accurate forces are required, PREC=Accurate can be used in combination with an increased energy cutoff.
                continue

            #	PRECFOCK
            #	PROUTINE
            #	PSTRESS
            #	PSUBSYS
            #	PTHRESHOLD
            #	QMAXFOCKAE
            #	QSPIRAL
            #	QUAD EFG
            #	RANDOM SEED
            #	ROPT
            #	RWIGS
            #	SAXIS
            #	SCSRAD
            #	SHAKEMAXITER
            #	SHAKETOL

            # SIGMA - specifies the width of the smearing (eV)
            if key=="SIGMA":
                dictionary[key] = float( value )
                continue

            #	SMASS
            #	SMEARINGS
            #	SPRING
            #	STEP MAX
            #	STEP SIZE
            #	SYMPREC
            # SYSTEM - just a user comment

            if key=="SYSTEM":
                # The "title string" defined by SYSTEM is for the user only.
                # Helps user identify what to do with the input file.
                continue

            #	TEBEG
            #	TEEND
            #	TIME
            #	TSUBSYS
            #	VALUE MAX
            #	VALUE MIN
            #	VCUTOFF
            #	VDW A1
            #	VDW A2
            #	VDW C6
            #	VDW CNRADIUS
            #	VDW D
            #	VDW R0
            #	VDW RADIUS
            #	VDW S6
            #	VDW S8
            #	VDW SR
            #	VOSKOWN
            #	WC
            #	WEIMIN
            #	ZVAL
            warning("processINCAR() boo boo with {0}.".format(key),"(ignoring key)")
    return dictionary







































# ################################# RUN CODE ################################# #


# This lets us place functions after they're called.
# Essentially everything gets loaded before runtime.
# Note! This must append the script.
if __name__ == '__main__':
    INCAR = parseIncar( '3scf' )

    getTypeIncar( INCAR )
    printBriefIncar( INCAR )
    #printFullIncar( INCAR )
