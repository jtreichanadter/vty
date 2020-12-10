#!/usr/bin/env python
"""vty.py - Run setup for 'vty' utilities inside the python environment.
Abstract:
    Use this script to startup the interactive python environment, with vty
    tools loaded. It effectively imports all of the vty utilities along with
    providing helpful documentation.
Example:
    Should be run in sh via command 'python -B -i ~/.scripts/vty/vty.py'.
"""
# -----------------------------------------------------------------------------#
# -----------------------------------------------------------------------------#
# -----------------------------------------------------------------------------#



def sh( *command ):
    # doc string
    """Let's you run shell commands within the python environment.
    [0] Abstract:
            Runs shell command (as a new process), returning a print statement.
            All processes are relative to where the scripts (or environment)
            were started from.
            NOTE: commands are fairly limited to avoid extreme error cases.
            You can view these by running sh() without arguments.
    [1] Dependencies:
            import subprocess (run,check_output,STDOUT)
            import os (devnull)
    [1] Input:
        string OR string list - shell command
    [2] Output:
        string list - what normally prints to terminal, split by lines
        NOTE: If command throws an error, return is None
    """
    # list of commands we allow
    shellList = ['cat','cp','curl','cut','df','du','echo','find','grep','head','ls','more','mv','printf','pwd','rev','tac','tail','top','wc']
    # check for input
    try: command = command[0]
    except:
        print("\033[0;1;38;5;220mAllowed commands:\n\033[0;38;5;220m{0}\033[0m".format(shellList))
        return None
    # prune input
    if   isinstance(command,list) and isinstance(command[0],str): pass
    elif isinstance(command,str) and len(command.strip())==0: print( "\033[0;38;5;208mERROR: cannot run sh() with empty string\033[0m");   return None
    elif isinstance(command,str): command = command.split()
    else: print( "\033[0;38;5;208mERROR: cannot run sh() with input type: {0}\033[0m".format(type(command)) );   return None
    # filter by commands
    #breakIndices = [0].append([ idx+1 for idx,cmd in enumerate(command[:-1]) if cmd=="|" ])
    #print(breakIndices)
    #if not any( not cmd in shellList for cmd in command):
    if command[0] in shellList:
        # import statements
        from subprocess import run          as _sh
        from subprocess import check_output as _shCheck
        from subprocess import STDOUT       as _STDOUT
        from os         import devnull      as _DNULL
        # check for errors, then return
        with open(_DNULL,'w') as _FNULL:
            if _sh( command , stdout=_FNULL , stderr=_STDOUT ).returncode == 0:
                return str( _shCheck(command) )[1:].strip('\"\'\\n ').split('\\n')
        # for errors return blank
        return None
    else:
        print( "\033[0;38;5;208mERROR: sh() function only accepts commands in the following list:\n{0}\033[0m".format(shellList) );
        return None
def shell(*command): return sh( *command )





def info( *infoInput ):
    """Prints out function documentation on vty tools.
    [0] Abstract:
        To keep track of your work, this function concisely prints out
        your current variables. Note it ignores all functions and most
        dunder variables
    [1] Dependencies:
    [2] Input:
    [3] Output:

    """
    # CASE 1: input
    try:
        infoInput = infoInput[0]
        if not callable(infoInput):
            print("\033[0;1;38;5;166mCould not find that function, print list...\033[0m")
        description = str(infoInput.__doc__).strip('\n').strip()
        if description=="None":
            print("\033[1;38;5;220m{0}():   \033[0;1;38;5;166mNothing specified, sorry friend :(\033[0m\n".format(str(infoInput).split()[1]))
        else:
            # doctor the print statement
            firstLine = description[0:description.find('\n')]
            description = description.replace(firstLine,"\033[3m"+firstLine+"\033[23m")
            description = description.replace("Abstract:\n","\033[4mAbstract\033[24m:\n")
            description = description.replace("Input:\n","\033[4mInput\033[24m:\n")
            description = description.replace("Output:\n","\033[4mOutput\033[24m:\n")
            # print out
            print("\033[1;38;5;220m{0}():  \033[0;38;5;227m{1}\033[0m\n".format(str(infoInput).split()[1],description))
    # CASE 2: no input
    except:
        # list of things to hide while printing the workspace environment
        ignoreDunderList = ['__name__','__doc__','__package__','__loader__','__spec__','__annotations__','__builtins__','__cached__','__status__','__author__','__email__','__copyright__','__license__','__module__']
        ignoreFcnList = ['printWelcome','shell','ws','WS']
        # begin big print-out
        for key in sorted(globals()):
            # begin printing logic
            if key not in ignoreDunderList and key not in ignoreFcnList:
                try:
                    if not callable(globals()[key]): continue
                except:
                    if not hasattr(globals()[key],'__call__'): continue
                
                description = str(globals()[key].__doc__).strip('\n').strip()
                if description=="None": print("\033[1;38;5;220m{0}():  \033[0;1;38;5;166mNothing specified, sorry friend :(\033[0m".format(key))
                else:
                    # print only the 1st statement
                    firstLine = description[0:description.find('\n')]
                    print("\033[1;38;5;220m{0}():  \033[0;38;5;227m{1}\033[0m".format(key,firstLine))
    # do not return any value
    return






def workspace():
    """Prints global variables within the python environment.
    [0] Abstract:
        To keep track of your work, this function concisely prints out
        your current variables. Note it ignores all functions and most
        dunder variables
    [1] Dependencies:
        built-in 'type()' function
    [2] Input:
    [3] Output:

    """

    # list of things to hide while printing the workspace environment
    ignoreDunderList = ['__name__','__doc__','__package__','__loader__','__spec__','__annotations__','__builtins__','__cached__','__status__','__author__','__email__','__copyright__','__license__','__module__']

    # prune workspace input [also remove all functions]
    # (2 methods of fcn checking used for more universal compatibility)
    try:    workSpaceList = [ key for key in sorted(globals())   if   key not in ignoreDunderList  and  not callable(globals()[key])   ]
    except: workSpaceList = [ key for key in sorted(globals())   if   key not in ignoreDunderList  and  not hasattr(globals()[key],'__call__')   ]
    if len(workSpaceList)==0: return

    # print formatting value
    keySpace = max([ len(str(key)) for key in workSpaceList ])

    # begin printing loop
    for key in workSpaceList:
        # print the first variable bit
        print("\033[0;1;38;5;220m{0} ".format(str(key).rjust(keySpace)),end='')
        value = globals()[key]
        # run through cases...
        typeString = str(type(value)).split("'")[1]
        # CASE: integer, float, complex, or string
        if any([ isinstance(value,dataType) for dataType in [int,float,complex,str] ]):
            print("\033[0;38;5;227m{0}\033[0;2;38;5;227m {1}\033[0m".format(str(value),typeString))
        # CASE: list (iterate through dimensions, assume uniform list type)
        elif isinstance(value,list):
            listDimensions = [ ]
            tmpVal = value
            while isinstance(tmpVal,list):
                listDimensions.append( len(tmpVal) )
                tmpVal = tmpVal[0]
            typeList = str(type(tmpVal)).split("'")[1]
            print("\033[0;38;5;227mlist [{0}]\033[0;2;38;5;227m {1}\033[0m".format(  "x".join( str(l) for l in listDimensions)  ,  typeList  ))
        # CASE: dictionary (compare with defined types)
        elif isinstance(value,dict):

            print("\033[0;38;5;227mdictionary")
        else:
            print("\033[0;2;38;5;227m {1}\033[0m".format(str(value),typeString))
    print("\033[0m",end='')
def ws(): workspace()
def WS(): workspace()







### DEFINE ESSENTIAL DICTIONARIES
def makeINCARdictionary():
    """Returns an empty INCAR dictionary.
    [0] Abstract:
        This returns a python dictionary specific to INCAR data, which is an
        essential tool for processing VASP INCAR files via 'vty'.
    [1] Input:
    [2] Output:
        incar dictionary - keys are INCAR flags, all values are empty
    """
    _INCAR_LIST = [ "_filePath","ADDGRID" , "AEXX" , "AGGAC" , "AGGAX" , "ALDAC" , "ALGO" , "AMIN" , "AMIX" , "AMIX MAG" , "ANDERSEN PROB" , "ANTIRES" , "APACO" , "BMIX" , "BMIX MAG" , "CH LSPEC" , "CH NEDOS" , "CH SIGMA" , "CLL" , "CLN" , "CLNT" , "CLZ" , "CMBJ" , "CMBJA" , "CMBJB" , "CSHIFT" , "DEPER" , "DIMER DIST" , "DIPOL" , "DQ" , "EBREAK", "EDIFF", "EDIFFG", "EFIELD", "EFIELD PEAD" , "EINT" , "Electric Field Gradient" , "EMAX" , "EMIN" , "ENAUG" , "ENCUT" , "ENCUTFOCK" , "ENCUTGW" , "ENCUTGWSOFT" , "ENINI" , "EPSILON" , "EVENONLY" , "EVENONLYGW" , "FERDO" , "FERWE" , "FINDIFF" , "GGA" , "GGA COMPAT" , "HFLMAX" , "HFRCUT" , "HFSCREEN" , "HILLS BIN" , "HILLS H" , "HILLS W" , "HITOLER" , "I CONSTRAINED M" , "IALGO" ,"IBAND" , "IBRION" , "ICHARG" , "ICHIBARE" , "ICORELEVEL" , "IDIPOL" , "IEPSILON" , "IGPAR" , "IMAGES" , "IMIX" , "INCREM" , "INIMIX" , "INIWAV" , "IPEAD" , "ISIF" , "ISMEAR" , "ISPIN" , "ISTART" , "ISYM" , "IVDW" , "IWAVPR" , "KBLOCK" , "KGAMMA" , "KPAR" , "KPOINT BSE" , "KPUSE" , "KSPACING" , "LADDER" , "LAECHG" , "LAMBDA" , "LANGEVIN GAMMA" , "LANGEVIN GAMMA L" , "LASPH" , "LASYNC" , "LATTICE CONSTRAINTS" , "LBERRY" , "LBLUEOUT" , "LBONE" , "LCALCEPS" , "LCALCPOL" , "LCHARG" , "LCHIMAG" , "LCORR" , "LDAU" , "LDAUJ" , "LDAUL" , "LDAUPRINT" , "LDAUTYPE" , "LDAUU" , "LDIAG" , "LDIPOL" , "LEFG" , "LELF" , "LEPSILON" , "LFOCKAEDFT" , "LHARTREE" , "LHFCALC" , "LHYPERFINE" , "LKPROJ" , "LLRAUG" , "LMAXFOCK" , "LMAXFOCKAE" , "LMAXMIX" , "LMAXPAW" , "LMAXTAU" , "LMIXTAU" , "LMONO" , "LNABLA" , "LNMR SYM RED" , "LNONCOLLINEAR" , "LOCPROJ" , "LOPTICS" , "LORBIT" , "LORBMOM" , "LPARD" , "LPEAD" , "LPLANE" , "LREAL" , "LRPA" , "LSCAAWARE" , "LSCALAPACK" , "LSCALU" , "LSCSGRAD" , "LSELFENERGY" , "LSEPB" , "LSEPK" , "LSORBIT" , "LSPECTRAL" , "LSPECTRALGW" , "LSPIRAL" , "LSUBROT" , "LTHOMAS" , "LUSE VDW" , "LVDW EWALD" , "LVDW ONECELL" , "LVDWEXPANSION" , "LVHAR" , "LVTOT" , "LWANNIER90" , "LWANNIER90 RUN" , "LWAVE" , "LWRITE MMN AMN" , "LWRITE UNK" , "LWRITE WANPROJ" , "LZEROZ" , "M CONSTR" , "MAGMOM" , "MAXMEM" , "MAXMIX" , "MDALGO" , "METAGGA" , "MINROT" , "MIXPRE" , "ML FF AFILT2 MB" , "ML FF CDOUB" , "ML FF CSF" , "ML FF CSIG" , "ML FF CSLOPE" , "ML FF CTIFOR" , "ML FF EATOM" , "ML FF IAFILT2 MB" , "ML FF IBROAD1 MB" , "ML FF IBROAD2 MB" , "ML FF ICOUPLE MB" , "ML FF ICUT1 MB" , "ML FF ICUT2 MB" , "ML FF IERR" , "ML FF IREG MB" , "ML FF ISAMPLE" , "ML FF ISCALE TOTEN MB" , "ML FF ISOAP1 MB" , "ML FF ISOAP2 MB" , "ML FF ISTART" , "ML FF IWEIGHT" , "ML FF LAFILT2 MB" , "ML FF LBASIS DISCARD" , "ML FF LCONF DISCARD" , "ML FF LCOUPLE MB" , "ML FF LCRITERIA" , "ML FF LEATOM MB" , "ML FF LHEAT MB" , "ML FF LMAX2 MB" , "ML FF LMLFF" , "ML FF LMLMB" , "ML FF LNORM1 MB" , "ML FF LNORM2 MB" , "ML FF MB MB" , "ML FF MCONF" , "ML FF MCONF NEW" , "ML FF MHIS" , "ML FF MRB1 MB" , "ML FF MRB2 MB" , "ML FF MSPL1 MB" , "ML FF MSPL2 MB" , "ML FF NATOM COUPLED MB" , "ML FF NDIM SCALAPACK" , "ML FF NHYP1 MB" , "ML FF NHYP2 MB" , "ML FF NMDINT" , "ML FF NR1 MB" , "ML FF NR2 MB" , "ML FF NWRITE" , "ML FF RCOUPLE MB" , "ML FF RCUT1 MB" , "ML FF RCUT2 MB" , "ML FF SIGV0 MB" , "ML FF SIGW0 MB" , "ML FF SION1 MB" , "ML FF SION2 MB" , "ML FF W1 MB" , "ML FF W2 MB" , "ML FF WTIFOR" , "ML FF WTOTEN" , "ML FF WTSIF" , "NBANDS" , "NBANDSGW" , "NBANDSO" , "NBANDSV" , "NBLK" , "NBLOCK" , "NBMOD" , "NBSEEIG" , "NCORE" , "NCRPA BANDS" , "NDAV" , "NEDOS" , "NELECT" , "NELM" , "NELMDL" , "NELMIN" , "NFREE" , "NGX" , "NGXF" , "NGY" , "NGYF" , "NGYROMAG" , "NGZ" , "NGZF" , "NKRED" , "NKREDX" , "NKREDY" , "NKREDZ" , "NLSPLINE" , "NMAXFOCKAE" , "LMAXFOCKAE" , "NOMEGA" , "NOMEGAPAR" , "NOMEGAR" , "NPACO" , "NPAR" , "NPPSTR" , "NSIM" , "NSUBSYS" , "NSW" , "NTARGET STATES" , "NTAUPAR" , "NUPDOWN" , "NWRITE" , "ODDONLY" , "ODDONLYGW" , "OFIELD A" , "OFIELD KAPPA" , "OFIELD Q6 FAR" , "OFIELD Q6 NEAR" , "OMEGAMAX" , "OMEGAMIN" , "OMEGATL" , "PARAM1" , "PARAM2" , "PFLAT" , "PHON LBOSE" , "PHON LMC" , "PHON NSTRUCT" , "PHON NTLIST" , "PHON TLIST" , "PLEVEL" , "PMASS" , "POMASS" , "POTIM" , "PREC" , "PRECFOCK" , "PROUTINE" , "PSTRESS" , "PSUBSYS" , "PTHRESHOLD" , "QMAXFOCKAE" , "QSPIRAL" , "QUAD EFG" , "RANDOM SEED" , "ROPT" , "RWIGS" , "SAXIS" , "SCSRAD" , "SHAKEMAXITER" , "SHAKETOL" , "SIGMA" , "SMASS" , "SMEARINGS" , "SPRING" , "STEP MAX" , "STEP SIZE" , "SYMPREC" , "SYSTEM" , "TEBEG" , "TEEND" , "TIME" , "TSUBSYS" , "VALUE MAX" , "VALUE MIN" , "VCUTOFF" , "VDW A1" , "VDW A2" , "VDW C6" , "VDW CNRADIUS" , "VDW D" , "VDW R0" , "VDW RADIUS" , "VDW S6" , "VDW S8" , "VDW SR" , "VOSKOWN" , "WC" , "WEIMIN" , "ZVAL" ]
    return { key:None for key in _INCAR_LIST }
def makePOSCARdictionary():
    """Returns an empty POSCAR dictionary.
    [0] Abstract:
        This returns a python dictionary specific to POSCAR data, which is an
        essential tool for processing VASP POSCAR and CONTCAR files via 'vty'.
    [1] Input:
    [2] Output:
        poscar(contcar) dictionary - all values are empty
    """
    _POSCAR_LIST = [ "_filePath" , "A" , "B" , "comment" , "elements" , \
                     "positions" , "positionsDirect" , "selectiveDynamics" , \
                     "velocities" ]
    return { key:None for key in _POSCAR_LIST }
def makeKPOINTSdictionary():
    """Returns an empty KPOINTS dictionary.
    [0] Abstract:
        This returns a python dictionary specific to KPOINTS data, which is an
        essential tool for processing VASP KPOINTS files via 'vty'.
    [1] Input:
    [2] Output:
        kpoints dictionary - all values are empty
    """
    _KPOINTS_LIST = [ "_filePath" , "auto" , "comment" , "format" , "k" , \
                      "kx" , "ky" , "kz" , "path" , "style" , "weights" ]
    return { key:None for key in _KPOINTS_LIST }
def makeOUTCARdictionary():
    """Returns an empty OUTCAR dictionary.
    [0] Abstract:
        This returns a python dictionary specific to OUTCAR data, which is an
        essential tool for storing VASP OUTCAR file data via 'vty'.
    [1] Input:
    [2] Output:
        outcar dictionary - keys are various outcar data, all values are empty
    """
    _OUTCAR_LIST = [ "_filePath" , "bands" , "bandsFermi" , "bandNum" ,  \
                     "energyGS" , "fermi" , "forces" , "ionicSteps" , \
                     "ionNum" , "iterations" , "k" , "kNum" , "muX" , "muY" , \
                     "muZ" , "occupancies" , "positions" , "stepNum" , \
                     "timing" , "weights" ]
    return { key:None for key in _OUTCAR_LIST }
def dictionaryMatch( dictionary1 , dictionary2 ):
    """Compare fields of 2 python dictionaries.
    [0] Abstract:
        Return true only if both input dictionaries have identical fields.
    [1] Dependencies:
    [2] Input:
        dictionary - (generic)
        dictionary - (generic)
    [3] Output:
        boolean
    """
    return \
        not any( not key in dictionary1.keys() for key in dictionary2.keys() ) \
        and \
        not any( not key in dictionary2.keys() for key in dictionary1.keys() )



# IMPORT VTY TOOLS
from vtySH_findFile  import *













# -----------------------------------------------------------------------------#
# -----------------------------------------------------------------------------#
# -----------------------------------------------------------------------------#

__module__    = "VTY"
__status__    = "Development" # (Prototype,Development,Production)
__author__    = "(Jonathan) Tyler Reichanadter"
__email__     = "jtreichanadter@berkeley.edu"
__copyright__ = "Copyright 2020, Neaton Group at UC Berkeley"
__license__   = "GPL"





# --- SHELL SCRIPT ----------------------------------------------------------- #
if __name__ == '__main__':
    # temporarily assign ascii escape sequences
    reset  = "\033[0m"
    gold   = "\033[0;38;5;227m"
    goldi  = "\033[0;3;38;5;220m"
    GOLD   = "\033[0;1;38;5;220m"
    red    = "\033[0;38;5;208m"
    RED    = "\033[0;1;38;5;208m"
    # big welcome message
    print("{0}  WELCOME!{1}".format(GOLD,gold))
    print("You can now use {0}vty{1} utilities interactively!".format(GOLD,gold))
    print("--> Run {0}workspace(){1} or {0}ws(){1} to view currently loaded variables.".format(GOLD,gold))
    print("--> Run {0}sh(){1} via an input string to execute shell commands.".format(GOLD,gold))
    print("--> Run {0}info(){1} to list all available {0}vty{1} utilities.".format(GOLD,gold))
    print("--> {0}info(){1} also accepts script names and VASP lists for more details!".format(GOLD,gold))
    print("--> {0}exit(){1} takes you out of python, wiping all memory...".format(GOLD,gold))
    print("{0}Have fun playing in the sandbox :){1}".format(gold,reset))
    print("{0}\nNOTE:{1} This python environment has been launched in\n      {0}{ourDirectory}".format(goldi,gold,ourDirectory=sh("pwd")[0]))
    # check modules
    try:
        # try import statments
        from subprocess import run as _sh
        from os import devnull as _DNULL
        del _sh
        del _DNULL
    except:
        # print warning
        print("\n{RED}WARNING! CANNOT IMPORT subprocess and os FROM subprocess LIBRARY.{red}\nThis means we CANNOT run shell commands in this environment; BIG problem :({reset}".format(RED=RED,red=red,reset=reset))
        print("{0}Tragically this means that most {1}vty{0} utilities may not work...{2}\n".format(red,RED,reset))
    # close out
    print(reset)
    del reset,gold,goldi,GOLD,red,RED
