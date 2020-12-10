#! /usr/bin/env python





# ############################## MAIN EXECUTION ############################## #

def parseOutcar( *Directory ):

    # specify file-name from Directory
    try:
        Directory = Directory[0]
        fileName = generateFileName(Directory)
    except:
        fileName = generateFileName('.')

    # creates dictionary for storing 'outcar' file details
    outcar = { key:[] for key in keyListOutcar() }
    # parse staging variable
    stage = 0
    mu = "muX"

    ##### #  START LOOP # #####
    # parsing sequence reads down-to-up (required by 'stage' variable)
    with open(fileName,'r') as outcarFile:
        for line in outcarFile:

            # skip empty lines and separate words
            if len(line.strip())==0: continue
            words = line.split()

            # [0] normal operation (check flags for cases)
            if stage==0:
                # get fermi energy
                if len(words)>2 and words[0]=="E-fermi":
                    outcar["fermi"] = float( words[2] )
                # get ground-state energy
                elif "free energy" in line:
                    outcar["energyGS"].append( float(words[4]) )
                # get iteration number
                elif len(words)>3 and words[1]=="Iteration":
                    outcar["iterations"].append( int(float(words[3][:-1])) )
                # CASE 1: read k-points direct
                elif outcar["k"]==[] and words[0]=="k-points" and len(words)>5 and words[2]=="reciprocal":
                    stage=1
                    continue
                # CASE 2: read band structure
                elif len(words)>1 and words[0]=="band" and words[1]=="No.":
                    outcar["bands"].append([])
                    outcar["occupancies"].append([])
                    stage=2
                    continue
                # CASE 3: read magnetization
                elif len(words)==2 and words[0]=="magnetization":
                    if words[1]=="(x)": mu = "muX"
                    if words[1]=="(y)": mu = "muY"
                    if words[1]=="(z)": mu = "muZ"
                    outcar[mu].append([])
                    stage=3
                    continue
                # CASE 4: read positions and forces
                elif len(words)==3 and words[0]=="POSITION":
                    outcar["positions"].append([])
                    outcar["forces"].append([])
                    stage = 4
                    continue
                # CASE N: read final time info
                elif len(words)>6 and words[0]=="General" and words[1]=="timing":
                    stage = -1
                    continue

            # [1] reading out k-points
            if stage==1:
                if isNumberString(line) and len(words)>3:
                    outcar["k"].append([ float(words[0]) , float(words[1]) , float(words[2]) ])
                    outcar["weights"].append( float(words[3]) )
                elif not outcar["k"]==[]:
                    stage = 0
                continue

            # [2] read band energies
            if stage==2:
                if isNumberString(line) and len(words)>2:
                    outcar["bands"][-1].append(float(words[1]))
                    outcar["occupancies"][-1].append(float(words[2]))
                else:
                    stage = 0

            # [3] read magnetization
            if stage==3:
                if isNumberString(line) and len(words)>4:
                    outcar[mu][-1].append(float(words[4]))
                elif not outcar[mu][-1]==[]:
                    stage = 0

            # [4] read position
            if stage==4:
                if isNumberString(line) and len(words)==6:
                    tmpPositions = []
                    tmpForces = []
                    for j in range(3):
                        tmpPositions.append(float(words[j]))
                        tmpForces.append(float(words[j+3]))
                    outcar["positions"][-1].append(tmpPositions)
                    outcar["forces"][-1].append(tmpForces)
                elif not outcar["positions"][-1]==[]:
                    stage = 0

            # [N] Closing statement
            if stage==-1:
                if words[0]=="Elapsed" and words[1]=="time":
                    outcar["timing"] = float( words[3] )


    # process and return the 'outcar' dictionary
    outcar = processOutcar( outcar )
    return outcar






# ############################## PRINTING DATA ############################## #


# publish data to new files
def publishOutcar( dictionary , *Directory ):

    # check that we are looking at a OUTCAR dictionary
    assert isDictionaryOutcar( dictionary ), "publishOutcar() input was not an outcar dictionary!"

    # set writing dictionary (default to current directory)
    try:
        Directory = Directory[0]
        Directory = Directory.strip('/')
    except:
        Directory = "."

    # data references
    BANDS = dictionary["bands"]
    EF = dictionary["fermi"]
    OCC = dictionary["occupancies"]
    K = dictionary["k"]
    ITER = dictionary["iterations"]
    EGS = dictionary["energyGS"]
    POS = dictionary["positions"]
    FORCE = dictionary["forces"]
    # range dimensions
    I = range(len(K))
    II = range(len(BANDS))
    J = range(len(BANDS[0]))
    L = range(len(ITER))
    M = range(dictionary["ionNum"])
    if len(K)<len(BANDS): N = int( dictionary["stepNum"] )

    # make scaled k-index
    kIDX = [ 0.0 for i in I ]
    for i in I[1:]: kIDX[i] = kIDX[i-1] + ( (K[i][0]-K[i-1][0])**2 + (K[i][1]-K[i-1][1])**2 + (K[i][2]-K[i-1][2])**2  )**(1/2)
    kIDX = [ kIDX[i]/kIDX[-1] for i in I ]


    # [1] bands-Ef VS kpts
    #
    # col(0) = k-index
    # col(1) = k vector, X-component
    # col(2) = k vector, Y-component
    # col(3) = k vector, Z-component
    # col(4,N+4) = N bands at k-index
    #
    # try statement succeeds when it's a relaxation calculation,
    # makes a band-structure file for each ionic step
    try:
        for n in range(N):
            # WRITE DATA
            fileName = Directory + "/.OUTCAR_BANDS_s" + str(n) + ".dat"
            open(fileName,"w").write("".join( "{0} {1} {2} {3} {4} {5}\n".format(i+1,kIDX[i],K[i][0],K[i][1],K[i][2]," ".join(str(BANDS[i+n*len(K)][j]-EF) for j in J) ) for i in I ));
    except:
        # WRITE DATA
        fileName = Directory + "/.OUTCAR_BANDS.dat"
        open(fileName,"w").write("".join( "{0} {1} {2} {3} {4} {5}\n".format(i+1,kIDX[i],K[i][0],K[i][1],K[i][2]," ".join(str(BANDS[i][j]-EF) for j in J) ) for i in I ))
    # WRITE META FILE
    fileName = Directory + "/.OUTCAR_BANDS_info.txt"
    open(fileName,"w").write("This file contains each band's energy (eV) with respect to the Fermi energy at each k-point\n")
    open(fileName,"a").write("Total of "+str(len(K))+" k-points (rows) and "+str(len(BANDS[0]))+" bands")
    open(fileName,"a").write("\nCOLUMN 001 :   k-point index")
    open(fileName,"a").write("\nCOLUMN 002 :   scaled k-point indices by k-vector magnitude, normalized [0,1]")
    open(fileName,"a").write("\nCOLUMN 003 :   k-point X-component")
    open(fileName,"a").write("\nCOLUMB 004 :   k-point Y-component")
    open(fileName,"a").write("\nCOLUMB 005 :   k-point Z-component")
    open(fileName,"a").write("".join("\nCOLUMN {0} :   energy of band number {1}".format(str(j+6).zfill(3),j+1) for j in J))


    # [2] occupancies VS kpts
    #   col(0) = k-index
    #   col(1) = normalized k-index
    #   col(1) = k vector, X-component
    #   col(2) = k vector, Y-component
    #   col(3) = k vector, Z-component
    #   col(4,J+4) = J bands at a k-index
    try:
        for n in range(N):
            # WRITE DATA
            fileName = Directory + "/.OUTCAR_OCCUPANCIESs" + str(n+1) + ".dat"
            open(fileName,"w").write("".join( "{0} {1} {2} {3} {4} {5}\n".format(i+1,kIDX[i],K[i][0],K[i][1],K[i][2]," ".join(str(OCC[i+n*len(K)][j]) for j in J) ) for i in I ));
    except:
        # WRITE DATA
        fileName = Directory + "/.OUTCAR_OCCUPANCIES.dat"
        open(fileName,"w").write("".join( "{0} {1} {2} {3} {4} {5}\n".format(i+1,kIDX[i],K[i][0],K[i][1],K[i][2]," ".join(str(OCC[i][j]) for j in J) ) for i in I ))
    # WRITE META FILE
    fileName = Directory + "/.OUTCAR_OCCUPANCIES_info.txt"
    open(fileName,"w").write("This file contains each band's occupancy at each k-point\n")
    open(fileName,"a").write("Total of "+str(len(K))+" k-points (rows) and "+str(len(BANDS[0]))+" bands")
    open(fileName,"a").write("\nCOLUMN 001 :   k-point index")
    open(fileName,"a").write("\nCOLUMN 002 :   scaled k-point indices by k-vector magnitude, normalized to range [0,1]")
    open(fileName,"a").write("\nCOLUMN 003 :   k-point X-component")
    open(fileName,"a").write("\nCOLUMB 004 :   k-point Y-component")
    open(fileName,"a").write("\nCOLUMB 005 :   k-point Z-component")
    open(fileName,"a").write("".join("\nCOLUMN {0} :   occupancy of band number {1}".format(str(j+6).zfill(3),j+1) for j in J))



    # [3] GS-energy vs iteration
    #
    # col(0) = ionic step
    # col(1) = iteration
    # col(2) = ground-state energy
    STEP = [ 1 for l in L ]
    for l in L[1:]:
        if ITER[l]<ITER[l-1]: STEP[l:] = [ step+1 for step in STEP[l:] ]
    # WRITE DATA
    fileName = Directory + "/.OUTCAR_GROUNDSTATES.dat"
    open(fileName,"w").write("".join( "{0} {1} {2}\n".format(STEP[l],ITER[l],EGS[l]) for l in L ))
    # WRITE META FILE
    fileName = Directory + "/.OUTCAR_GROUNDSTATES_info.txt"
    open(fileName,"w").write("This file contains the electronic ground-state energy after each iteration\n")
    open(fileName,"a").write("Total of "+str(len(L))+" iterations and "+str(STEP[-1]-1)+" ionic steps")
    open(fileName,"a").write("\nCOLUMN 1 :   ionic step index")
    open(fileName,"a").write("\nCOLUMN 2 :   electronic step index")
    open(fileName,"a").write("\nCOLUMB 3 :   ground-state energy (eV)")


    # POSTIONS vs steps
    #
    # col(0) = ionic step
    # col(1) = ground-state energy
    # col(2,M+1) = position-x [over ions]
    # col(M+2,2M+1) = position-y [over ions]
    # col(2M+2,3M+1) = position-z [over ions]
    #
    # col(0) = ionic step
    # col(1) = ground-state energy
    # col(2,M+1) = force-x [over ions]
    # col(M+2,2M+1) = force-y [over ions]
    # col(2M+2,3M+1) = force-z [over ions]
    # col(3M+2,4M+1) = force-mag [over ions]
    #
    # only for relaxation calculations
    try:
        # trips the try statement
        N = N+0;
        # make gsEnergy steps
        GS = []
        for l in L[1:]:
            if ITER[l]<ITER[l-1]: GS.append(EGS[l-1])
        GS.append(EGS[-1])
        # WRITE DATA: POSITIONS
        fileName = Directory + "/.OUTCAR_POSITIONS.dat"
        open(fileName,"w").write("".join( "{0} {1} {2}\n".format( \
                n+1,str(GS[n]) , \
                " ".join(str(POS[n][m][0])+" "+str(POS[n][m][1])+" "+str(POS[n][m][2]) for m in M) ) \
                for n in range(N) ));
        # WRITE META: POSITIONS
        fileName = Directory + "/.OUTCAR_POSITIONS_info.txt"
        open(fileName,"w").write("This file contains position information (in cartesian coordinates, Angstroms) after each ionic step\n")
        open(fileName,"a").write("Total of "+str(N)+" ionic steps")
        open(fileName,"a").write("\nCOLUMN 001 :   ionic step index")
        open(fileName,"a").write("\nCOLUMN 002 :   ground-state energy (eV)")
        open(fileName,"a").write("".join("\nCOLUMN {1} :   X-position of ion {0}\nCOLUMN {2} :   Y-position of ion {0}\nCOLUMN {3} :   Z-position of ion {0}".format(m+1,str(3*m+3).zfill(3),str(3*m+4).zfill(3),str(3*m+5).zfill(3)) for m in M ))
        # make force norm array
        FORCE_T = [ [ 0.0 for m in M ] for n in range(N) ]
        for n in range(N):
            for m in M:
                FORCE_T[n][m] = ( FORCE[n][m][0]**2 + FORCE[n][m][1]**2 + FORCE[n][m][2]**2 )**(1/2)
        # WRITE DATA: FORCES
        fileName = Directory + "/.OUTCAR_FORCES.dat"
        open(fileName,"w").write("".join( "{0} {1} {2}\n".format( \
                n+1,str(GS[n]) , \
                " ".join(str(FORCE[n][m][0])+" "+str(FORCE[n][m][1])+" "+str(FORCE[n][m][2])+" "+str(FORCE_T[n][m]) for m in M) ) \
                for n in range(N) ));
        # WRITE META: FORCES
        fileName = Directory + "/.OUTCAR_FORCES_info.txt"
        open(fileName,"w").write("This file contains force information (in standard units, Ev/Angstrom) after each ionic step\n")
        open(fileName,"a").write("Total of "+str(N)+" ionic steps")
        open(fileName,"a").write("\nCOLUMN 001 :   ionic step index")
        open(fileName,"a").write("\nCOLUMN 002 :   ground-state energy (eV)")
        open(fileName,"a").write("".join("\nCOLUMN {1} :   X-force of ion {0}\nCOLUMN {2} :   Y-force of ion {0}\nCOLUMN {3} :   Z-force of ion {0}\nCOLUMN {4} :   |force| of ion {0}".format( m+1 , str(4*m+3).zfill(3) , str(4*m+4).zfill(3) , str(4*m+5).zfill(3) , str(4*m+6).zfill(3) ) for m in M ))

    except: pass














# ############################# HELPER-FUNCTIONS ############################# #



# acquire fileName
def generateFileName(Directory):
    # strip off tail
    Directory = Directory.rstrip('/')
    # base file name
    NAME = "OUTCAR"
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
        # try other names
        for j in range(10):
            try:
                checkFile = open(fileName+"-"+str(j),'r')
                checkFile.close()
                warning( "Could not open file '{0}'".format(fileName) , "Using {0}".format(fileName+"-"+str(j)) )
                return fileName+"-"+str(j)
            except:
                pass
    # could not find the file
    warning("Failed to find {0} file in parse{0}.py!!".format(NAME))
    return NAME



# checks if the given dictionary is generated from an INCAR file
def isDictionaryOutcar( dictionary ):
    if not isinstance(dictionary,dict): return False
    return \
        not any( not key in keyListOutcar() for key in dictionary.keys() ) \
        and \
        not any( not key in dictionary.keys() for key in keyListOutcar() )


# prints a full OUTCAR dictionary vertically
def printFullOutcar( dictionary ):
    # check that we are looking at a OUTCAR dictionary
    assert isDictionaryOutcar( dictionary )
    # prints everything (including empty keys)
    for item in dictionary.items():
        print("{0}:  {1}".format(item[0],item[1])  )


# prints the nonempty elements of a OUTCAR dictionary vertically
def printOutcar( dictionary ):
    # check that we are looking at a OUTCAR dictionary
    assert isDictionaryOutcar( dictionary )
    # prints everything (including empty keys)
    for item in dictionary.items():
        if not item[1]==None:
            if isinstance(item[1],list):
                print("{0}:".format(item[0]))
                if len(item[1])>20:
                    printMat(item[1][0:20])
                    print("."); print("."); print(".")
                else:
                    printMat(item[1])
            else:
                print("{0}: {1}".format(item[0],item[1]))

# prints the element types of an OUTCAR dictionary
def printBriefOutcar( dictionary ):
    # check that we are looking at a OUTCAR dictionary
    assert isDictionaryOutcar( dictionary )
    # prints everything (including empty keys)
    for item in dictionary.items():
        if item[1]==None:
            print( "{0}:  EMPTY".format(item[0]) )
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


# verifies that the provided string is composed purely of numbers
def isNumberString( fullString ):
    # iterate through words and catch any non-float values
    for word in fullString.split():
        try: tmp = float(word)
        except: return False
    # return true if we get here
    return True


# prints OUTCAR parser warnings
def warning(*argv):
    print("\033[1;33m--- WARNING : parseOutcar ---\033[0;33m")
    for arg in argv:
        print(arg)
    print("\033[0m")


# prints out matrices
def printMat( matrix ):
    # case: scalar:
    if not isinstance(matrix,list):
        print(matrix)
        return
    # case: row vector
    if len(matrix)==1:
        print(matrix[0])
        return
    # otherwise:
    for i in range(len(matrix)):
        # case: column vector
        if not isinstance(matrix[0],list):
            print("|",matrix[i],"|")
        # case: matrix
        else: print(matrix[i])









# ############################# OUTCAR-HELPERS ############################## #


# generate list of OUTCAR dictionary keys
def keyListOutcar():

    return [ "bands" , "bandsFermi" , "bandNum" ,  "energyGS" , "fermi" , "forces" , \
             "ionicSteps" , "ionNum" , "iterations" , "k" , "kNum" , \
             "muX" , "muY" , "muZ" , "occupancies" , "positions" , \
             "stepNum" , "timing" , "weights" ]


# process the raw OUTCAR info
def processOutcar( dictionary ):

    # reduce magnetic moment components to their last lines (if calculation is magnetic)
    try:
        dictionary["muX"] = dictionary["muX"][-1]
        dictionary["muY"] = dictionary["muY"][-1]
        dictionary["muZ"] = dictionary["muZ"][-1]
    except:
        pass

    # fill in a couple of missed dictionary keys
    dictionary["kNum"] = len(dictionary["k"])
    dictionary["bandNum"] = len(dictionary["bands"][0])
    dictionary["stepNum"] = len(dictionary["positions"])
    dictionary["ionNum"] = len(dictionary["positions"][0])

    # create ionic step number array
    for i in range(1,len(dictionary["iterations"])):
        if dictionary["iterations"][i]<=dictionary["iterations"][i-1]:
            dictionary["ionicSteps"].append(dictionary["iterations"][i-1])
    dictionary["ionicSteps"].append(dictionary["iterations"][-1])

    # create fermi-zeroed bandstructure
    I = range(dictionary["kNum"])
    J = range(dictionary["bandNum"])
    dictionary["bandsFermi"] = [[ dictionary["bands"][i][j] - dictionary["fermi"] for j in J ] for i in I]

    # reset any empty entries
    for item in dictionary.items():
        if item[1]==[]: dictionary[item[0]] = None

    # return processed 'outcar'
    return dictionary









# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

__module__    = "VTY"
__status__    = "Development" # (Prototype,Development,Production)
__author__    = "(Jonathan) Tyler Reichanadter"
__email__     = "jtreichanadter@berkeley.edu"
__copyright__ = "Copyright 2020, Neaton Group at UC Berkeley"
__license__   = "GPL"




# ################################# RUN CODE ################################# #

if __name__ == '__main__':

    OUTCAR = parseOutcar( '.' )
    printBriefOutcar(OUTCAR)
    #publishOutcar(OUTCAR)
    print(OUTCAR["muZ"])
