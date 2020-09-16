#! /usr/bin/env python



# ############################## MAIN EXECUTION ############################## #
def parsePoscar( *Directory ):

    # specify file-name from Directory
    try:
        Directory = Directory[0]
        fileName = generateFileName(Directory)
    except:
        fileName = generateFileName('.')

    # creates dictionary for storing 'incar' file details
    poscar = { key:None for key in poscarKeyList() }


    with open(fileName,'r') as poscarFile:

        # create staging index, for knowing what to parse
        stage = 0
        # index for assigning array values
        idx = 0

    	##### #  START LOOP # #####
        # parsing sequence reads down-to-up (required by 'stage' variable)
        for line in poscarFile:

            # depricate comments, skip empty lines
            if stage>0:
                try: line = line[ 0:line.index("#") ]
                except ValueError: pass
                if len(line.strip())==0: continue
            line = line.strip()

            # (6) assign volicity values
            if stage==6:
                vector = line.split()
                # avoid inappropriate input line
                if isNumberString(line) and len(vector)==3:
                    for j in range(3): V[idx][j] = float(vector[j])
                    idx = idx + 1
                if idx==total:
                    poscar["velocities"] = V
                    idx = 0

            # (5) assign position values
            if stage==5:
                vector = line.split()
                # avoid inappropriate input line
                if isNumberString(line) and len(vector)==3:
                    for j in range(3): P[idx][j] = float(vector[j])
                    idx = idx + 1
                if idx==total:
                    if position_cartesian:
                        poscar["positions"] = P
                    else:
                        poscar["positionsDirect"] = P
                    idx = 0
                    # initialize velocity array
                    V = [ [ None for j in range(3) ] for i in range(total) ]
                    # step forward
                    stage = 6


            # (4) identify lattice format (and check for selective dynamics)
            if stage==4:
                words = line.split()
                character = words[0][0]
                if character in ['s','S']:
                    selectiveDynamics = True
                    continue
                elif character in ['c','C','k','K']:
                    position_cartesian = True
                else:
                    position_cartesian = False
                poscar["selectiveDynamics"] = selectiveDynamics
                # step forward
                stage = 5


            # (3) form ion list
            if stage==3:
                # catch ion labels (if provided)
                if not isNumberString(line):
                    labels = line.split()
                    continue
                # count ions
                vector = line.split()
                total = 0
                for term in vector:
                    total = total + int(float( term ))
                # label ions (if appropriate)
                try:
                    if len(labels)==len(vector):
                        elements = []
                        for i in range(len(vector)):
                            for j in range(int(float(vector[i]))):
                                elements.append(labels[i])
                        poscar["elements"] = elements
                # EXIT
                finally:
                    # initialize position array (and selective dynamics bool)
                    P = [ [ None for j in range(3) ] for i in range(total) ]
                    selectiveDynamics = False
                    # step forward
                    stage = 4

            # (2) Get lattice vectors
            if stage==2:
                # verify we are looking at lattice vectors
                vector = line.split()
                if not isNumberString(line) or len(vector)!=3: continue
                # pull lattice components
                for j in range(3): A[idx][j] = scaling * float( vector[j] )
                idx = idx + 1
                # complete lattice assignment
                if idx==3:
                    poscar["A"] = A
                    idx = 0
                    # step forward
                    stage = 3

            # (1) Get scaling factor
            if stage==1:
                scaling = abs(float( line ))
                # initialize lattice matrix
                A = [ [ None for j in range(3) ] for i in range(3) ]
                # step forward
                stage = 2

            # (0) Read comment header
            if stage==0:
                poscar["comment"] = line.strip()
                # step forward
                stage = 1



    # post-process raw POSCAR dictionary after parsing
    poscar = processPoscar(poscar)

    # DONE
    return poscar






















# ############################# HELPER-FUNCTIONS ############################# #




# define pi
pi = float("3.1415926535897932384626433832795")



# prints POSCAR parser warnings
def warning(*argv):
    print("\033[1;33m--- WARNING : parsePOSCAR.py ---\033[0;33m")
    for arg in argv:
        print(arg)
    print("\033[0m")


# acquire fileName
def generateFileName(Directory):
    # strip off tail
    Directory = Directory.rstrip('/')
    # base file name
    NAME = "POSCAR"
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



# checks if the given dictionary is generated from a POSCAR file
def isDictionaryPoscar( dictionary ):
    return \
        not any( not key in poscarKeyList() for key in dictionary.keys() ) \
        and \
        not any( not key in dictionary.keys() for key in poscarKeyList() )


# verifies that the provided string is composed purely of numbers
def isNumberString( fullString ):
    # iterate through words and catch any non-float values
    for word in fullString.split():
        try: tmp = float(word)
        except: return False
    # return true if we get here
    return True


# prints the nonempty elements of an POSCAR dictionary vertically
def printPoscar( dictionary ):
    # check that we are looking at a POSCAR dictionary
    assert isDictionaryPoscar( dictionary ) , "'printPoscar()' READ AN INAPPROPRIATE DICTIONARY"
    # prints everything (including empty keys)
    for item in dictionary.items():
        if not item[1]==None:
            if isinstance(item[1],list):
                print("{0}:".format(item[0]))
                printMat(item[1])
            else:
                print("{0}: {1}".format(item[0],item[1]))


# prints a full POSCAR dictionary vertically
def printFullPoscar( dictionary ):
    # check that we are looking at an POSCAR dictionary
    assert isDictionaryPoscar( dictionary ) , "'printFullPoscar()' READ AN INAPPROPRIATE DICTIONARY"
    # prints everything (including empty keys)
    for item in dictionary.items():
            if isinstance(item[1],list):
                print("\033[1;37m{0}:\033[0m".format(item[0]))
                printMat(item[1])
            else:
                print("\033[1;37m{0}:\033[0m {1}".format(item[0],item[1]))

# prints the element types of an POSCAR dictionary
def printBriefPoscar( dictionary ):
    # check that we are looking at a POSCAR dictionary
    assert isDictionaryPoscar( dictionary )
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


# return list of periodic elements
def elementList():
    return [ "H"  , "He" , \
             "Li" , "Be" , \
             "B"  , "C"  , "N"  , "O"  , "F"  , "Ne" , \
             "Na" , "Mg" , \
             "Al" , "Si" , "P"  , "S"  , "Cl" , "Ar" , \
             "K"  , "Ca" , \
             "Sc" , "Ti" , "V"  , "Cr" , "Mn" , "Fe" , "Co" , "Ni" , "Cu" , "Zn" , \
             "Ga" , "Ge" , "As" , "Se" , "Br" , "Kr" , \
             "Rb" , "Sr" , \
             "Y"  , "Zr" , "Nb" , "Mo" , "Tc" , "Ru" , "Rh" , "Pd" , "Ag" , "Cd" , \
             "In" , "Sn" , "Sb" , "Te" , "I"  , "Xe" , \
             "Cs" , "Ba" , \
             "La" , "Ce" , "Pr" , "Nd" , "Pm" , "Sm" , "Eu" , "Gd" , "Tb" , "Dy" , "Ho" , "Er" , "Tm" , "Yb" , \
             "Lu" , "Hf" , "Ta" , "W"  , "Re" , "Os" , "Ir" , "Pt" , "Au" , "Hg" , \
             "Tl" , "Pb" , "Bi" , "Po" , "At" , "Rn" , \
             "Fr" , "Ra" , \
             "Ac" , "Th" , "Pa" , "U"  , "Np" , "Pu" , "Am" , "Cm" , "Bk" , "Cf" , "Es" , "Fm" , "Md" , "No" , \
             "Lr" , "Rf" , "Db" , "Sg" , "Bh" , "Hs" , "Mt" , "Ds" , "Rg" , "Cn" , \
             "Nh" , "Fl" , "Mc" , "Lv" , "Ts" , "Og" ]


# check rank of input (level of list nesting, ignores dictionary,tuple,)
def checkRank( input ):
    # toss out non-lists (base case)
    if not isinstance(input,list): return 0
    # recursively tally the rank
    if len(input)>1: return 1 + checkRank(input[0])
    # continue tally, but exclude length-1 dimensions!
    return checkRank(input[0])


# take scalar multiplication on vector
def scaMult( scalar , tensor ):
    # handle case where arguments are flipped
    if isinstance(scalar,list) and not isinstance(tensor,list): scaMult(tensor,scalar)
    # check tensor rank
    rank = checkRank(tensor)
    assert rank<3, "'scaMult()' TENSOR RANK EXCEEDS LIMIT!"
    # case: matrix
    if rank==2:
        result = [ [ 0.0 for j in range(len(tensor[0])) ] for i in range(len(tensor)) ]
        for i in range(len(tensor)):
            for j in range(len(tensor[0])):
                result[i][j] = scalar * tensor[i][j]
        return result
    # case: row vector
    try:
        result = [ [ 0.0 for j in range(tensor[0]) ] for i in range(1) ]
        for j in range(len(tensor[0])): result[0][j] = scalar*tensor[0][j]
        return result
    # case: column vector
    except:
        result = [ 0.0 for i in range(len(tensor)) ]
        for i in range(len(tensor)): result[i] = scalar*tensor[i]
        return result


# take dot product
def vecDot( vec1 , vec2 ):
    # check for appropriate vectors
    assert type(vec1)==type(vec2), "'vecDot()' TRYING TO USE OFF-TYPE VECTORS!"
    assert len(vec1)==len(vec2), "'vecDot()' TRYING TO USE OFF-LENGTH VECTORS!"
    # compute dot product
    dotSum = 0
    for j in range(len(vec1)): dotSum = dotSum + vec1[j]*vec2[j]
    return dotSum


# take cross product
def vecCross( vec1 , vec2 ):
    # check for appropriate vectors
    assert type(vec1)==type(vec2), "'vecCross()' TRYING TO USE OFF-TYPE VECTORS!"
    assert len(vec1)==len(vec2), "'vecCross()' TRYING TO USE OFF-LENGTH VECTORS!"
    assert len(vec1)==3, "'vecCross()' TRYING TO USE NON-3D VECTORS!"
    # compute cross product
    cross = [ 0 for j in range(3) ]
    cross[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1]
    cross[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2]
    cross[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0]
    return cross


# form vector projection (v1,v2) v1
def vecProj( vec1 , vec2 ):
    # check for appropriate vectors
    assert type(vec1)==type(vec2), "'vecProj()' TRYING TO USE OFF-TYPE VECTORS!"
    assert len(vec1)==len(vec2), "'vecProj()' TRYING TO USE OFF-LENGTH VECTORS!"
    # compute projection
    scaling = vecDot(vec1,vec2) / vecDot(vec1,vec1)
    projection = [ 0 for j in range(len(vec1)) ]
    for j in range(len(vec1)): projection[j] = scaling * vec1[j]
    return projection

# take transpose (accepts matrices, vectors, and scalars)
def transpose( matrix ):
    # limit tensor rank
    rank = checkRank(matrix)
    assert rank<3 , "'transpose()' TENSOR RANK EXCEEDS 2"
    # case: matrix
    if rank==2:
        # initialize matrix
        newMat = [ [ 0.0 for j in range(len(matrix)) ] for i in range(len(matrix[0])) ]
        # populate matrix
        for i in range(len(matrix)):
            for j in range(len(matrix[0])):
                newMat[j][i] = matrix[i][j]
        return newMat
    # case: vector
    elif rank==1:
        # case: nested vector
        try:
            # initialize vector
            newVec = [ 0 for i in range(len(matrix[0])) ]
            # populate vector
            for j in range(len(matrix[0])):
                newVec[j] = matrix[0][j]
            return newVec
        # case: simple vector
        except:
            # initialize vector
            newVec = [ [0 for j in range(len(matrix)) ] for i in range(1) ]
            # populate vector
            for i in range(len(matrix)):
                newVec[0][i] = matrix[i]
            return newVec
    # case: scalar
    else:
        return matrix

# take matrix product (can involve vectors)
def product( mat1 , mat2 ):
    # verify no higher-order tensors
    assert checkRank(mat1)<3 , "'product()' RANK OF INPUT1 EXCEEDS LIMIT"
    assert checkRank(mat2)<3 , "'product()' RANK OF INPUT2 EXCEEDS LIMIT"
    # get dimensions
    try: M = len(mat1)
    except: M = 1
    try: K = len(mat1[0])
    except: K = 1
    try: K2 = len(mat2)
    except: K2 = 1
    try: N = len(mat2[0])
    except: N = 1
    # handle scalar exceptions
    if (M==1 and K==1): return product(mat2,mat1)
    # check for compatible dimensions
    assert K==K2 , "'product()' INCOMPATIBLE DIMENSIONS"
    # BEGIN CASES
    if N>1:
        # case: (MxK) (KxN) --> (MxN) matrix
        if M>1:
            newMat = [ [ 0 for j in range(N) ] for i in range(M) ]
            for i in range(M):
                for j in range(N):
                    if K>1:
                        for k in range(K):
                            newMat[i][j] = newMat[i][j] + mat1[i][k] * mat2[k][j]
                    else: newMat[i][j] = newMat[i][j] + mat1[i] * mat2[0][j]
        # case: (1xK) (KxN) --> (1xN) row vector
        else:
            newMat = [ [ 0 for j in range(N) ] for i in range(1) ]
            for j in range(N):
                if K>1:
                    for k in range(K):
                        newMat[0][j] = newMat[0][j] + mat1[0][k] * mat2[k][j]
                else: newMat[0][j] = newMat[0][j] + mat1[0] * mat2[0][j]
    else:
        # case: (MxK) (Kx1) --> (Mx1) column vector
        if M>1:
            newMat = [ 0 for i in range(M) ]
            for i in range(M):
                if K>1:
                    for k in range(K):
                        newMat[i] = newMat[i] + mat1[i][k] * mat2[k]
                else: newMat[i] = newMat[i] + mat1[i] * mat2
        # case (1xK) (Kx1) --> (1x1) scalar
        else:
            if K>1:
                newMat = 0
                for k in range(K):
                    newMat = newMat + mat1[0][k] * mat2[k]
            else: newMat = mat1 * mat2
    # RETURN
    return newMat

# prints out matrices
def printMat( matrix ):
    # verify no higher-order tensors
    assert checkRank(matrix)<3 , "'printMat()' RANK EXCEEDS LIMIT"
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






# ############################## POSCAR-HELPERS ############################## #



# returns an alphabetical list of poscar keys
def poscarKeyList():
    return [ "A" , "B" , "comment" , "elements" , "positions" , \
             "positionsDirect" , "selectiveDynamics" , "velocities" ]


# updates raw POSCAR dictionary (after parsing)
def processPoscar( dictionary ):

    # check that we are looking at an POSCAR dictionary
    assert isDictionaryPoscar( dictionary ) , "'processPoscar()' READ AN INAPPROPRIATE DICTIONARY"

    # pull lattice (for computations)
    A = dictionary["A"]

    # (B) - compute reciprocal lattice vectors
    B = [ [ 0 for i in range(3) ] for j in range(3) ]
    B[0] = vecCross( A[1] , A[2] )
    B[0] = scaMult(2*pi/vecDot(A[0],B[0]),B[0])
    B[1] = vecCross( A[2] , A[0] )
    B[1] = scaMult(2*pi/vecDot(A[1],B[1]),B[1])
    B[2] = vecCross( A[0] , A[1] )
    B[2] = scaMult(2*pi/vecDot(A[2],B[2]),B[2])
    dictionary["B"] = B

    # (positions) - fill both cartesian and direct coordinate forms
    if dictionary["positions"] != None:
        # cartesian --> direct
        dictionary["positionsDirect"] = product(dictionary["positions"],scaMult(1/2/pi,transpose(B)))
    elif dictionary["positionsDirect"] != None:
        # direct --> cartesian
        dictionary["positions"] = product(dictionary["positionsDirect"],A)
    else:
        print("ERROR: NO POSITION DATA COLLECTED?")

    # (elements) - check element input
    elements = dictionary["elements"]
    if elements == None:
        dictionary["elements"] = [ "??" for i in range(len(dictionary["positions"])) ]
    else:
        for i in range(len(elements)):
            if elements[i] in elementList(): dictionary["elements"][i] = elements[i].rjust(2)

    # (velocities) - fill in empty item with zero components
    if dictionary["velocities"] == None:
        dictionary["velocities"] = [ [ 0.0 for j in range(3) ] for i in range(len(dictionary["positions"])) ]

    # DONE
    return dictionary




















# ################################# RUN CODE ################################# #


# This lets us place functions after they're called.
# Essentially everything gets loaded before runtime.
# Note! This must append the script.
if __name__ == '__main__':
    POSCAR = parsePoscar( '2relax/' )
    printBriefPoscar( POSCAR )
