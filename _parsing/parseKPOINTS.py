#! /usr/bin/env python





# ############################## MAIN EXECUTION ############################## #
def parseKpoints( *Directory ):

    # specify file-name from Directory
    try:
        Directory = Directory[0]
        fileName = generateFileName(Directory)
    except:
        fileName = generateFileName('.')

    # creates dictionary for storing 'kpoints' file details
    kpoints = { key:None for key in kpointsKeyList() }
    # assign initial kpoint styling (may change due to text flags)
    kpoints["style"] = "mesh"
    kpoints["format"] = "direct"
    kset = []
    weights = []
    # variables for extra k-point file controls
    readTetragonal = False
    shiftMesh = False
    # create staging index, for knowing what to parse
    stage = 0


    ##### #  START LOOP # #####
    # parsing sequence reads down-to-up (required by 'stage' variable)
    with open(fileName,'r') as kpointsFile:
        for line in kpointsFile:


            # depricate comments, skip empty lines
            if stage>0:
                try: line = line[ 0:line.index("!") ]
                except ValueError: pass
                if len(line.strip())==0: continue
                line = line.strip()

            # [3] read in the path or mesh
            if stage==3:
                # execute only if they're numbers
                if isNumberString(line):
                    # --- urgent control lines --- #
                    if readTetragonal: # Reading tetragonal connections
                        # STUB
                        warning("'readTetragonal' has not yet been coded!")
                        pass
                    elif shiftMesh: # Shifting to auto-grid
                        # STUB
                        warning("'shiftMesh' has not yet been coded!")
                        pass
                    # --- standard control --- #
                    else:
                        values = line.split()
                        # append new k-point
                        kset.append([0.0,0.0,0.0])
                        # populate new k-point
                        for j in range(3):
                            kset[len(kset)-1][j] = float(values[j])
                        # fill in the weight
                        try: weights.append(float(values[3]))
                        except: weights.append(1.0)
                else:
                    # reset urgent controls
                    readTetragonal = False
                    shiftMesh = False
                    # check for new controls ...
                    stage = 2

            # [2] Check format of k-points (skipped if auto)
            if stage==2:
                # read first character
                flag = line[0]
                # Is auto-type with Gamma?
                if   flag in ['G','g']: kpoints["auto"] = "gamma"
                # Is auto-type with Monkhorst Pack?
                elif flag in ['M','m']: kpoints["auto"] = "monkhorstPack"
                # Is path-type?
                elif flag in ['L','l']: kpoints["style"] = "line"
                # Is cartesian?
                elif flag in ['C','c','K','k']: kpoints["format"] = "cartesian"
                # Check for tetrahedral reader
                elif flag in ['T','t']: readTetragonal = True
                # Check for shifted auto-grid
                elif flag in ['V','v']: shiftMesh = True
                # assume numeric lines follow...
                stage = 3

            # [1] grab k-point number specifier
            if stage==1:
                kpoints["path"] = int(float( line ))
                # check if auto-style is specified
                if kpoints["path"] == 0:
                    kpoints["style"] = "auto"
                    kpoints["auto"] = "default"
                # step forward
                stage = 3

            # [0] Read header (always a comment)
            if stage==0:
                kpoints["comment"] = line.strip()
                # step forward
                stage = 1


    # fill in the 'k' set of the dictionary
    kpoints["k"] = kset
    kpoints["weights"] = weights

    # process and return the 'kpoints' dictionary
    kpoints = processKPOINTS( kpoints )
    return kpoints









# ############################# HELPER-FUNCTIONS ############################# #


# acquire fileName
def generateFileName(Directory):
    # strip off tail
    Directory = Directory.rstrip('/')
    # base file name
    NAME = "KPOINTS"
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
def isDictionaryKpoints( dictionary ):
    return \
        not any( not key in kpointsKeyList() for key in dictionary.keys() ) \
        and \
        not any( not key in dictionary.keys() for key in kpointsKeyList() )


# prints a full KPOINTS dictionary vertically
def printFullKpoints( dictionary ):
    # check that we are looking at a KPOINTS dictionary
    assert isDictionaryKpoints( dictionary )
    # prints everything (including empty keys)
    for item in dictionary.items():
        print("{0}:  {1}".format(item[0],item[1])  )


# prints the nonempty elements of a KPOINTS dictionary vertically
def printKpoints( dictionary ):
    # check that we are looking at a KPOINTS dictionary
    assert isDictionaryKpoints( dictionary )
    # prints everything (including empty keys)
    for item in dictionary.items():
        if not item[1]==None:
            if isinstance(item[1],list):
                print("{0}:".format(item[0]))
                printMat(item[1])
            else:
                print("{0}: {1}".format(item[0],item[1]))

# prints the element types of an POSCAR dictionary
def printBriefKpoints( dictionary ):
    # check that we are looking at a POSCAR dictionary
    assert isDictionaryKpoints( dictionary )
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


# prints KPOINTS parser warnings
def warning(*argv):
    print("\033[1;33m--- WARNING : parseKpoints ---\033[0;33m")
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









# ############################# KPOINTS-HELPERS ############################## #



# generate list of KPOINTS dictionary keys
def kpointsKeyList():
    # keys:
    #   auto    -   special flag for 'automatic' kpoint mesh formatting
    #                 (either 'default', 'gamma', or 'monkhorstPack')
    #                 (nonempty only when 'style'=='auto')
    #   comment -   header comment
    #   format  -   either 'cartesian' or 'direct' coordinates for 'k'
    #                 ('cartesian' will need the reciprocal lattice from POSCAR)
    #   k       -   list of k-points in VASP default format
    #   kx      -   3D array k mesh for the x-dimension
    #   ky      -   3D array k mesh for the y-dimension
    #   kz      -   3D array k mesh for the z-dimension
    #   path    -   number of k-points between paths, or total point number
    #   style   -   How to read points 'mesh', 'auto', or 'line'
    #   weights -   How each k-point is weighted in calculations
    return [ "auto" , "comment" , "format" , "k" , "kx" , "ky" , "kz" , "path" , "style" , "weights" ]


# process the raw KPOINTS info
def processKPOINTS( dictionary ):

    # --- process 'k' data ---
    # [A] Manual (default) Input
    if dictionary["style"]=="mesh":
        warning("file 'parseKpoints.py' hasn't been coded for","the 'mesh' setting")
        pass
    # [B] Automatic Generation
    elif dictionary["style"]=="auto":
        # [B1] auto: defaut-mode
        if dictionary["auto"]=="default":
            warning("file 'parseKpoints.py' hasn't been coded for","the 'full auto' setting")
            pass
        # [B2] auto: gamma-centered mode
        elif dictionary["auto"]=="gamma":
            Dx = int(dictionary["k"][0][0])
            Dy = int(dictionary["k"][0][1])
            Dz = int(dictionary["k"][0][2])
            k = [ [0.0,0.0,0.0] for j in range(Dx*Dy*Dz) ]
            for s in range(Dx*Dy*Dz):
                u = s % Dx
                v = (s // Dx) % Dy
                w = (s // (Dx*Dy)) % Dz
                k[s][0] = u/Dx
                if k[s][0]>0.5: k[s][0] = k[s][0] - 1.0
                k[s][1] = v/Dy
                if k[s][1]>0.5: k[s][1] = k[s][1] - 1.0
                k[s][2] = w/Dz
                if k[s][2]>0.5: k[s][2] = k[s][2] - 1.0
            dictionary["k"] = k
        # [B3] auto: Monkhorst Pack mode
        elif dictionary["auto"]=="monkhorstPack":
            warning("file 'parseKpoints.py' hasn't been coded for","the 'auto monkhorstPack' setting")
            pass
    # [C] Path Input
    elif dictionary["style"]=="line":
        # prune k-paths
        kset = dictionary["k"]
        del kset[1:len(kset)-1:2]
        # forge new k-points array
        pts = dictionary["path"]
        paths = len(kset)-1
        k = [ [0.0,0.0,0.0] for j in range(paths*pts) ]
        # populate kpoints array
        for path in range(paths):
            for i in range(pts):
                for j in range(3):
                    k[path*pts+i][j] = ((pts-1-i)/(pts-1)) * kset[path][j] + ((i)/(pts-1)) *kset[path+1][j]
        dictionary["k"] = k


    # populate 'k,x,y,z' meshes (if relevant)
    if not dictionary["style"]=="line":
        pass


    # normalize weights
    if not dictionary["style"]=="mesh":
        dictionary["weights"] = None
    else:
        sum = 0.0
        weights = dictionary["weights"]
        for weight in weights:
            sum = sum + weight
        for j in range(len(weights)):
            weights[j] = weights[j]/sum


    # return processed 'kpoints'
    return dictionary











# ################################# RUN CODE ################################# #


# This lets us place functions after they're called.
# Essentially everything gets loaded before runtime.
# Note! This must append the script.
if __name__ == '__main__':
    KPOINTS = parseKpoints( '4band/KPOINTS' )
    print(" ")
    printKpoints(KPOINTS)
