#! /usr/bin/env python


# ############################## MAIN EXECUTION ############################## #
def parseProcar( *Directory ):

    # specify file-name from Directory
    try:
        Directory = Directory[0]
        fileName = generateFileName(Directory)
    except:
        fileName = generateFileName('.')

    # creates dictionary for storing 'procar' file details
    procar = { key:[] for key in keyListProcar() }
    # parsing variable
    stage = 0
    q = 0


    ##### #  START LOOP # #####
    # parsing sequence reads down-to-up (required by 'stage' variable)
    with open(fileName,'r') as procarFile:
        for line in procarFile:

            # skip empty lines and separate words
            if len(line.strip())==0 and stage!=3: continue
            words = line.split()

            # [3] check for the complex projection scheme
            if stage==3:
                if len(line.strip())==0:
                    # check exiting condition
                    if j<procar["bandNum"]-1:
                        stage = 1
                        continue
                    # exit back to k-points
                    else:
                        stage = 0
                        continue
                if len(words)>0 and words[0]=="ion":
                    k = 0
                elif isNumberString(line) and len(words)>=10:
                    if int(float(words[0]))==k:
                        # assign imaginary part
                        k = int(float(words[0]))
                        for s in range(1,10): procar["projectionsC"][i][j][k-1][s-1] = complex( procar["projectionsC"][i][j][k-1][s-1].real , float(words[s]) )
                    else:
                        # assign real part
                        k = int(float(words[0]))
                        for s in range(1,10): procar["projectionsC"][i][j][k-1][s-1] = complex( float(words[s]) , 0.0 )
                else: stage = 0

            # [0] normal operation (check for cases)
            if stage==0:
                # check for new k-point
                if words[0]=="k-point":
                    # increment k index
                    i = int(float( words[1] )) - 1
                    # USUAL OPERATION
                    if len(words)==9:
                        # append new k-point
                        procar["k"][i] = [ float(words[3]) , float(words[4]) , float(words[5]) ]
                        procar["weights"][i] = float(words[8])
                    # WHEN THE K-VECTOR STICKS (due to - signs)
                    else:
                        # acquire k-vector string
                        vector = "".join( (" "+word) for word in words[3:-3])
                        vector = vector.strip()
                        # identify vector components via nonnumeric indices
                        term = 'x'
                        idx = 0
                        for l in range(1,len(vector)):
                            if term=='x' and not vector[l].isdigit() and vector[l]!='.':
                                xComponent = float(vector[:l-1].strip())
                                idx = l
                                term='y'
                            elif term=='y' and l>idx and not vector[l].isdigit() and vector[l]!='.':
                                yComponent = float(vector[idx:l-1].strip())
                                idx = l
                                term='z'
                        zComponent = float(vector[idx:].strip())
                        # assign values
                        procar["k"][i] = [ xComponent , yComponent , zComponent ]
                        procar["weights"][i] = float(words[-1])
                    # move to next stage
                    stage = 1
                # assign meta variables (should occur first, and just once)
                elif len(words)==12 and words[0]=="#":
                    # assign variables
                    procar["kNum"] = int(float(words[3]))
                    procar["bandNum"] = int(float(words[7]))
                    procar["ionNum"] = int(float(words[11]))
                    # initialize arrays
                    procar["weights"] = [ -1.0 for i in range(procar["kNum"]) ]
                    procar["k"] = [ [ -1.0 for j in range(3) ] for i in range(procar["kNum"]) ]
                    procar["bands"] = [ [ -1.0 for j in range(procar["bandNum"]) ] for i in range(procar["kNum"])]
                    procar["occupancies"] = [ [ -1.0 for j in range(procar["bandNum"]) ] for i in range(procar["kNum"])]
                    procar["projectionsTotal"] = [ [ [ -1.0 for k in range(10) ] for j in range(procar["bandNum"]) ] for i in range(procar["kNum"]) ]
                    procar["muProjectionXTotal"] = [ [ [ -1.0 for k in range(10) ] for j in range(procar["bandNum"]) ] for i in range(procar["kNum"]) ]
                    procar["muProjectionYTotal"] = [ [ [ -1.0 for k in range(10) ] for j in range(procar["bandNum"]) ] for i in range(procar["kNum"]) ]
                    procar["muProjectionZTotal"] = [ [ [ -1.0 for k in range(10) ] for j in range(procar["bandNum"]) ] for i in range(procar["kNum"]) ]
                    procar["projections"] = [ [ [ [ -1.0 for l in range(10) ] for k in range(procar["ionNum"]) ] for j in range(procar["bandNum"]) ] for i in range(procar["kNum"]) ]
                    procar["projectionsC"] = [ [ [ [ complex(-0.0,-0.0) for l in range(9) ] for k in range(procar["ionNum"]) ] for j in range(procar["bandNum"]) ] for i in range(procar["kNum"]) ]
                    procar["muProjectionX"] = [ [ [ [ -1.0 for l in range(10) ] for k in range(procar["ionNum"]) ] for j in range(procar["bandNum"]) ] for i in range(procar["kNum"]) ]
                    procar["muProjectionY"] = [ [ [ [ -1.0 for l in range(10) ] for k in range(procar["ionNum"]) ] for j in range(procar["bandNum"]) ] for i in range(procar["kNum"]) ]
                    procar["muProjectionZ"] = [ [ [ [ -1.0 for l in range(10) ] for k in range(procar["ionNum"]) ] for j in range(procar["bandNum"]) ] for i in range(procar["kNum"]) ]


            # [1] entered a k-point
            if stage==1:
                # staring point for a band reading
                if len(words)>7 and words[0]=="band":
                    # increment band index
                    j = int(float( words[1] )) - 1
                    # append to bandstructure and occupancies
                    procar["bands"][i][j] = float(words[4])
                    procar["occupancies"][i][j] = float(words[7])
                    # move to projection stage
                    stage = 2


            # [2] entered a band-idx
            if stage==2:
                # append new values
                if isNumberString(line) and len(words)>10:
                    # increment ion index
                    k = int(float( words[0] )) - 1
                    # populate projection
                    #procar[proj][i][j][k] = [ float(words[1]) , float(words[2]) , float(words[3]) , float(words[4]) , float(words[5]) , float(words[6]) , float(words[7]) , float(words[8]) , float(words[9]) , float(words[10]) ]
                    procar[proj][i][j][k] = [ float(words[s]) for s in range(1,11)  ]
                # check for start
                elif words[0]=="ion":
                    proj = "projections"
                    continue
                # check magnetic element
                elif words[0]=="tot":
                    # populate projection total
                    #procar[proj+"Total"][i][j] = [ float(words[1]) , float(words[2]) , float(words[3]) , float(words[4]) , float(words[5]) , float(words[6]) , float(words[7]) , float(words[8]) , float(words[9]) , float(words[10]) ]
                    procar[proj+"Total"][i][j] = [ float(words[s]) for s in range(1,11)  ]
                    # identify new stuff
                    if proj=="projections": proj = "muProjectionX"
                    elif proj=="muProjectionX": proj = "muProjectionY"
                    elif proj=="muProjectionY": proj = "muProjectionZ"
                    # finished reading
                    else: stage = 3

    # return empty dictionary values to 'None' type
    for element in procar.items():
        if element[1]==[]: procar[element[0]] = None

    # process and return the 'procar' dictionary
    procar = processProcar( procar )
    return procar


# ############################## PRINTING DATA ############################## #
def publishProcar( dictionary , *Directory ):

    # check that we are looking at a OUTCAR dictionary
    assert isDictionaryProcar( dictionary ), "\n\npublishProcar() input was not an outcar dictionary!\n"

    # set writing dictionary (default to current directory)
    try:
        Directory = Directory[0]
        Directory = Directory.rstrip('/')
    except:
        Directory = "."

    # data references
    K = dictionary["k"]
    BANDS = dictionary["bands"]
    OCC = dictionary["occupancies"]
    # range dimensions
    I = range(len(K))
    J = range(len(BANDS[0]))
    M = range(dictionary["ionNum"])
    # make scaled k-index
    kIDX = [ 0.0 for i in I ]
    for i in I[1:]: kIDX[i] = kIDX[i-1] + ( (K[i][0]-K[i-1][0])**2 + (K[i][1]-K[i-1][1])**2 + (K[i][2]-K[i-1][2])**2  )**(1/2)
    kIDX = [ kIDX[i]/kIDX[-1] for i in I ]

#    [ "bands" , "bandNum" , "bandsFermiIndex" , "ionNum" , "k" , \
#      "kNum" , "muProjectionX" , "muProjectionXTotal" , "muProjectionY" , \
#      "muProjectionYTotal" , "muProjectionZ" , "muProjectionZTotal" , \
#      "occupancies" , "projectionLabels" , "projections" , "projectionsTotal" , "weights" ]

    # [1] bands-Ef VS kpts
    #   col(0) = k-index
    #   col(1) = normalized k-index
    #   col(2) = k vector, X-component
    #   col(3) = k vector, Y-component
    #   col(4) = k vector, Z-component
    #   col(5,N+5) = N bands at k-index
    # WRITE DATA
    fileName = Directory + "/.PROCAR_BANDS.dat"
    open(fileName,"w").write("".join( "{0} {1} {2} {3} {4} {5}\n".format(i+1,kIDX[i],K[i][0],K[i][1],K[i][2]," ".join(str(BANDS[i][j]) for j in J) ) for i in I ))
    # WRITE META FILE
    fileName = Directory + "/.PROCAR_BANDS_info.txt"
    file = open(fileName,"w")
    file.write("This file contains each band's energy (eV) with respect to the Fermi energy at each k-point\n")
    file.write("Total of "+str(len(K))+" k-points (rows) and "+str(len(BANDS[0]))+" bands")
    file.write("\nCOLUMN 001 :   k-point index")
    file.write("\nCOLUMN 002 :   scaled k-point indices by k-vector magnitude, normalized [0,1]")
    file.write("\nCOLUMN 003 :   k-point X-component")
    file.write("\nCOLUMB 004 :   k-point Y-component")
    file.write("\nCOLUMB 005 :   k-point Z-component")
    file.write("".join("\nCOLUMN {0} :   energy of band number {1}".format(str(j+6).zfill(3),j+1) for j in J))



    # [2] muProjectionZTotal VS k-point/band
    MuZ = [ [ dictionary["muProjectionZTotal"][i][j][-1] for j in J ] for i in I ]
    # WRITE DATA
    fileName = Directory + "/.PROCAR_MUZ.dat"
    open(fileName,"w").write("".join( "{0} {1} {2} {3} {4} {5}\n".format(i+1,kIDX[i],K[i][0],K[i][1],K[i][2]," ".join(str(MuZ[i][j]) for j in J) ) for i in I ))



    # [3] ion projection max VS k-point/band
    # WRITE DATA
    for n in range(len( dictionary["projections"][0][0] )):
        fileName = Directory + "/.PROCAR_PROJECTION_TOTAL_ION" + str(n+1) + ".dat"
        open(fileName,"w").write("".join( "{0} {1} {2} {3} {4} {5}\n".format(i+1,kIDX[i],K[i][0],K[i][1],K[i][2]," ".join(str(dictionary["projections"][i][j][n][-1]) for j in J) ) for i in I ))









# ############################# HELPER-FUNCTIONS ############################# #


# acquire fileName
def generateFileName(Directory):
    # strip off tail
    Directory = Directory.rstrip('/')
    # base file name
    NAME = "PROCAR"
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
def isDictionaryProcar( dictionary ):
    return \
        not any( not key in keyListProcar() for key in dictionary.keys() ) \
        and \
        not any( not key in dictionary.keys() for key in keyListProcar() )


# prints a full PROCAR dictionary vertically
def printFullProcar( dictionary ):
    # check that we are looking at a PROCAR dictionary
    assert isDictionaryProcar( dictionary )
    # prints everything (including empty keys)
    for item in dictionary.items():
        print("{0}:  {1}".format(item[0],item[1])  )


# prints the nonempty elements of a PROCAR dictionary vertically
def printProcar( dictionary ):
    # check that we are looking at a PROCAR dictionary
    assert isDictionaryProcar( dictionary )
    # prints everything (including empty keys)
    for item in dictionary.items():
        if not item[1]==None:
            if isinstance(item[1],list):
                print("{0}:".format(item[0]))
                if len(item[1])>10:
                    printMat(item[1][0:10])
                    print("."); print("."); print(".")
                else:
                    printMat(item[1])
            else:
                print("{0}: {1}".format(item[0],item[1]))


# prints the element types of an PROCAR dictionary
def printBriefProcar( dictionary ):
    # check that we are looking at a PROCAR dictionary
    assert isDictionaryProcar( dictionary )
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


# prints PROCAR parser warnings
def warning(*argv):
    print("\033[1;33m--- WARNING : parseProcar ---\033[0;33m")
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









# ############################# PROCAR-HELPERS ############################## #


# generate list of PROCAR dictionary keys
def keyListProcar():
    # keys:
    #   bands           - band structure energy, 2D list.
    #                       ex. bands[3][16] == 4th kpt, 17th band
    #                       (idx 1 = Kpoint)
    #                       (idx 2 = Band)
    #   bandNum        - number of bands computed, float
    #   bandFermi       - index of fermi band (at each 'k'), 1D list.
    #                       (idx 1 = Kpoint)
    #   ionNum          - number of ions, float.
    #   k               - k-vectors, 2D list.
    #                       ex. k[4][1] == 5th kpt, y-component.
    #                       (idx 1 = Kpoint)
    #                       (idx 2 = vector component)
    #   kNum            - number of k-points, float.
    #   muProjectionX   - spin x-projection on wvfcns, 3D list.
    #                       ex. muX[89][63][2][3] == 90th kpt, 64th bnd, 3rd ion, 4th shell px
    #                       (idx 1 = Kpoint)
    #                       (idx 2 = Band)
    #                       (idx 3 = Ion)
    #                       (idx 4 = Orbital)
    #   muProjectionY   - the Y equivalent of 'muProjectionX'
    #   muProjectionZ   - the Z equivalent of 'muProjectionX'
    #   occupancies
    #   projections
    #   weights

    return [ "bands" , "bandNum" , "bandsFermiIndex" , "ionNum" , "fermi" , \
     "k" , "kNum" , "muProjectionX" , "muProjectionXTotal" , "muProjectionY" , \
     "muProjectionYTotal" , "muProjectionZ" , "muProjectionZTotal" , "occupancies" , \
     "projectionLabels" , "projections" , "projectionsC" , "projectionsTotal" , "weights" ]


# process the raw PROCAR info
def processProcar( dictionary ):

    # determine fermi index at each k-point.
    # (defined by the first occupancy below 0.5)
    OCC = dictionary["occupancies"]
    I = range(dictionary["kNum"])
    J = range(dictionary["bandNum"])
    # initialize
    dictionary["bandsFermiIndex"] = [ -1 for i in I ]
    # identify indices
    for i in I:
        for j in J:
            if dictionary["bandsFermiIndex"][i]==-1 and OCC[i][j]<0.5:
                dictionary["bandsFermiIndex"][i] = j

    # include projection labels
    dictionary["projectionLabels"] = [" s "," py"," pz"," px","dxy","dyz","dz2","dxz","dx2","tot"]

    # re-map complex
    proj  = dictionary["projections"]
    projC = dictionary["projectionsC"]
    dictionary["projectionsC"] = [[[[ complex( \
                    proj[q][b][i][s] * projC[q][b][i][s].real / (abs(projC[q][b][i][s])+1e-12) , \
                    proj[q][b][i][s] * projC[q][b][i][s].imag / (abs(projC[q][b][i][s])+1e-12) ) \
                for s in range(9) ] \
                for i in range(len(proj[0][0])) ] \
                for b in range(len(proj[0])) ] \
                for q in range(len(proj)) ]

    # return processed 'procar'
    return dictionary









# ################################# RUN CODE ################################# #
if __name__ == '__main__':
    PROCAR = parseProcar( '4band' )

    printBriefProcar(PROCAR)
    publishProcar(PROCAR)
