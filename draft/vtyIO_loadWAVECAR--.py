#!/usr/bin/env python
"""Extracts data from the WAVECAR Binary.
[0] Abstract:
    Parse data from
[1] Dependencies:

[2] Input:
    wavecar dictionary - specifying
[3] Output:
    wavecar dictionary
"""
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #




def loadWAVECAR( *wavecar ):

    # INPUT
    kpoint = 80
    band = 30
    # read file
    try: wavecar = wavecar[0]
    except:
        # no input
        print("Warning! Either there was no input into loadWAVECAR(), or \
        the input was not a valid WAVECAR dictionary")
        wavecar = { key:None for key in _keyListWavecar() }
        return wavecar




    # CASE 1 : string (filePath) input
    if isinstance(wavecar,str):
        with open(wavecar, mode='rb' ) as file:

            # --- 0 --- prepare content
            print("Extracting WAVECAR content...",end="")
            wavecarContent = file.read()
            from struct import unpack
            print("...content read.")

            # --- 1 --- read 1st header
            OUTPUT = unpack("ddd", wavecarContent[:24])
            print("\033[0;1;38;5;236mtmp = \033[0;38;5;239m{0} \033[2m({1})\033[0m".format(OUTPUT,type(OUTPUT[0])))
            recordLength = int(OUTPUT[0])
            spinIndex = int(OUTPUT[1])
            precision = int(OUTPUT[2])
            print("\033[0;38;5;089mrecordLength = {0}\nspinIndex = {1}\nprecision = {2}".format(recordLength,spinIndex,precision))

            # -- 1.5 -- identify next non-zero
            idx = 24
            OUTPUT = 0.0
            while OUTPUT==0.0:
                OUTPUT = float(unpack("d",wavecarContent[8*idx:8*idx+8])[0])
                idx += 1
            nextStep = idx-1
            print("\033[0m--> NEXT STEP: {0}".format(nextStep))

            # --- 2 --- read 2nd header
            OUTPUT = unpack( "dddddddddddd" , wavecarContent[ nextStep*8:nextStep*8+8*12 ] )
            print("\033[0;1;38;5;236mtmp = \033[0;38;5;239m{0} \033[2m({1})\033[0m".format(OUTPUT,type(OUTPUT[0])))
            kLength = int(OUTPUT[0])
            bandLength = int(OUTPUT[1])
            energyCutoff = int(OUTPUT[2])
            a1 = list(OUTPUT[3:6])
            a2 = list(OUTPUT[6:9])
            a3 = list(OUTPUT[9:12])
            print("\033[0;38;5;089mk = {0}\nb = {1}\nENCUT = {2}\n\033[0;38;5;134ma1 = {3}\na2 = {4}\na3 = {5}".format(kLength,bandLength,energyCutoff,a1,a2,a3))

            # --- 3 --- acquire reciprocal lattice
            ( b1 , b2 , b3 ) = vecReciprocals( a1 , a2 , a3 )
            print("\033[0;38;5;105mb1 = {0}\nb2 = {1}\nb3 = {2}\033[0m".format(b1,b2,b3))

            # --- 4 --- repeat wavetrans setup process
            from math import acos as arccos
            from math import asin as arcsin
            from math import sin
            from math import pi
            kMAX = (0.26246583099999998*energyCutoff)**0.5
            deg = 180/pi
            deg = 1
            # [A] plane-wave counting
            phi_12 =     arccos( vecDot(b1,b2)       / vecNorm(b1) / vecNorm(b2) )
            sinphi_123 = vecDot( vecCross(b1,b2),b3) / vecNorm(b3) / vecNorm(vecCross(b1,b2) )
            nb1maxA = int( 1  +  kMAX / vecNorm(b1) / abs(sin(phi_12)) )
            nb2maxA = int( 1  +  kMAX / vecNorm(b2) / abs(sin(phi_12)) )
            nb3maxA = int( 1  +  kMAX / vecNorm(b3) / abs( sinphi_123) )
            npmaxA = int( 4*pi/3 * nb1maxA*nb2maxA*nb3maxA + 0.5)
            print('\033[38;5;106mphi_12 = {0}\nsinphi_123 = {1}\nnpmaxA = {2}'.format( phi_12*deg , sinphi_123 , npmaxA ))
            # [B] plane-wave counting
            phi_13 =     arccos( vecDot(b1,b3)       / vecNorm(b1) / vecNorm(b3) )
            sinphi_123 = vecDot( vecCross(b1,b3),b2) / vecNorm(b2) / vecNorm(vecCross(b1,b3) )
            nb1maxB = int( 1  +  kMAX / vecNorm(b1) / abs(sin(phi_13)) )
            nb2maxB = int( 1  +  kMAX / vecNorm(b2) / abs( sinphi_123) )
            nb3maxB = int( 1  +  kMAX / vecNorm(b3) / abs(sin(phi_13)) )
            npmaxB = int( 4*pi/3 * nb1maxB*nb2maxB*nb3maxB + 0.5)
            print('\033[38;5;107mphi_13 = {0}\nsinphi_123 = {1}\nnpmaxB = {2}'.format( phi_13*deg , sinphi_123 , npmaxB ))
            # [C] plane-wave counting
            phi_23 =     arccos( vecDot(b2,b3)       / vecNorm(b2) / vecNorm(b3) )
            sinphi_123 = vecDot( vecCross(b2,b3),b1) / vecNorm(b1) / vecNorm(vecCross(b2,b3) )
            nb1maxC = int( 1  +  kMAX / vecNorm(b1) / abs( sinphi_123) )
            nb2maxC = int( 1  +  kMAX / vecNorm(b2) / abs(sin(phi_23)) )
            nb3maxC = int( 1  +  kMAX / vecNorm(b3) / abs(sin(phi_23)) )
            npmaxC = int( 4*pi/3 * nb1maxC*nb2maxC*nb3maxC + 0.5)
            print('\033[38;5;108mphi_23 = {0}\nsinphi_123 = {1}\nnpmaxC = {2}'.format( phi_23*deg , sinphi_123 , npmaxC ))
            # final G index limits
            nb1max = max([ nb1maxA , nb1maxB , nb1maxC ])
            nb2max = max([ nb2maxA , nb2maxB , nb2maxC ])
            nb3max = max([ nb3maxA , nb3maxB , nb3maxC ])
            print('\033[0;38;5;109mnb1max = {0}\nnb2max = {1}\nnb3max = {2}'.format( nb1max , nb2max , nb3max ))
            # 2x to handle two component spinors
            npmax = 2 * min([ npmaxA , npmaxB , npmaxC ])
            irec = 3 + (kpoint-1)*(bandLength+1)
            print('\033[0;38;5;110mnpmax = {0}\nirec = {1}\033[0m'.format( npmax , irec ))

            # --- 5 --- allocate memory
            igall = [ [ 0.0 for j in range(3) ] for i in range(npmax) ]
            coeff = [ 0.0 for i in range(npmax) ]



            # -- 5.5 -- identify first k-index
            # Done by tracking the 1st occupancy elements, and back-tracing
            idx = nextStep+12
            OUTPUT = 0.0
            while OUTPUT==0.0:
                OUTPUT = float(unpack("d",wavecarContent[8*idx:8*idx+8])[0])
                idx += 1
            nextStep = idx-1
            print("\033[0m--> NEXT STEP: {0}\n".format(nextStep))


            kIdx1 = 2
            kSpan = 65


            '''
            # -- 5.x -- testing stuff
            OUTPUT = [0.0,0.0,0.0]
            idx = nextStep
            q = 0
            Q = 0
            readingLimit = 0
            (vMIN,vMEAN,vMAX) = (1.0E5,0.0,-1.0E5)
            # BIG LOOP
            while readingLimit<1000000000:
                readingLimit+=1
                # pull data
                try: OUTPUT = list(unpack("ddddd",wavecarContent[idx*8:idx*8+5*8]))
                except: return
                print("idx-{0} ... ".format())
                OUTPUT = [float(OUTPUT[j]) for j in range(5)]
                # CASE 1: ending a value streak
                if OUTPUT[0]==0.0 and Q>0:
                    print("\033[38;5;226mvalues \033[2mx{0}".format(Q))
                    print("--> min  = {0}\n--> mean = {1}\n--> max  = {2}\033[0m".format(vMIN,vMEAN,vMAX))
                    (vMIN,vMEAN,vMAX) = (1.0E5,0.0,-1.0E5)
                    q=1; Q=0
                    idx += 1
                # CASE 2: continuing a zero streak
                elif OUTPUT[0]==0.0:
                    q+=1
                    idx += 1
                # CASE 3: k-point found
                elif isKpointListing(OUTPUT[:3]):
                    print("\033[38;5;160mkpoint = {0} \033[2mat {1}\033[0m".format(OUTPUT[:3],idx))
                    q=0; Q=0
                    idx += 3
                # CASE 4: some other header
                elif Q==0 and q>0 and any([ output==0.0 for output in OUTPUT ]) and any([ output!=0.0 for output in OUTPUT ]):
                    temp = OUTPUT[ : [ j for j in range(5) if OUTPUT[j]==0.0][0] ]
                    q=0; Q=0
                    temp2 = list(unpack("dddddddddd",wavecarContent[idx*8:idx*8+10*8]))[len(temp):]
                    temp2 = [ float(tempVal) for tempVal in temp2 ]
                    if any([ output!=0.0 for output in temp2 ]):
                        while len(temp)<5 and OUTPUT[len(temp)]==0.0:
                            temp.append(0.0)
                            q+=1
                    print("\033[38;5;226m{0}\033[0m".format(temp))
                    idx += len(temp)
                # CASE 5: ending zero streak
                elif OUTPUT[0]!=0.0 and q>0:
                    print("\033[0;38;5;248m0.0 \033[2mx{0}\033[0m".format(q))
                    Q=1; q=0
                    idx += 1
                # CASE 6: continuing a value streak
                elif OUTPUT[0]!=0.0:
                    Q+=1
                    if OUTPUT[0]<vMIN: vMIN = OUTPUT[0]
                    if OUTPUT[0]>vMAX: vMAX = OUTPUT[0]
                    vMEAN  =  ((Q-1)/Q) * vMEAN  +  (1/Q)*OUTPUT[0]
                    idx += 1
                # CASE 7: error...
                else:
                    print("ERROR: case not handled...")
                    print("q={0}  Q={1}  OUTPUT={2}".format(q,Q,OUTPUT))
                    return
            return
            '''



            return



            # --- 6 --- pull coefficients
            tmp = unpack( "dddd" , wavecarContent[ 8*irec:8*irec+4*8 ] )
            nplane = int(tmp[0])
            wk = tmp[1:]
            print("nplane = {0}\nwk = {1}".format(nplane,wk))





            return True

            print("\033[0;1;38;5;236mtmp = \033[0;38;5;239m{0}  \033[2m({1})\033[0m".format(tmp,type(tmp[0])))

            ncnt = 0
            sumkg = [ 0.0 , 0.0 , 0.0 ]
            for g3 in range(2*nb3max):
                g3p = g3
                if g3 > nb3max: g3p = g3 - 2*nb3max - 1
                for g2 in range(2*nb2max):
                    g2p = g2
                    if g2>nb2max: g2p = g2 - 2*nb2max - 1
                    for g1 in range(2*nb1max):
                        g1p = g1
                        if g1>nb1max: g1p = g1 - 2*nb1max - 1
                        for j in range(3):
                            sumkg[j] = (wk[0]+g1p) * b1[j]  +  (wk[1]+g2p) * b2[j]  +  (wk[2]+g3p) * b3[j]
                            gtot = vecDot(sumkg,sumkg) ** 0.5
                            etot = gtot**2 / 0.2624658310
                            if etot < energyCutoff:
                                ncnt += 1
                                igall[1][ncnt] = g1p
                                igall[2][ncnt] = g2p
                                igall[3][ncnt] = g3p





    # dumb return (for now)
    wavecar = { key:None for key in _keyListWavecar() }
    return wavecar



















# --- wavecar dictionary keys ------------------------------------------------ #

def _keyListWavecar():
    return [ "_filePath" , "A" , "B" , "b" , "chi" , "G" , "iSpin" , "k" , "precision" , "psi" , "recordLength" ]








# --- helper functions ------------------------------------------------------- #

def vecDot( vecA , vecB ):
    # check the inputs
    print( "\033[3;38;5;235m" , end='' )
    assert isinstance(vecA,list) and len(vecA)==3 and type(vecA[0]) in [int,float,complex], "\n\033[0;38;5;202mERROR: cross() input 1 must be a length-3 list of numbers, currently:\n{0}\033[0m".format(vecA)
    assert isinstance(vecB,list) and len(vecB)==3 and type(vecB[0]) in [int,float,complex], "\n\033[0;38;5;202mERROR: cross() input 2 must be a length-3 list of numbers, currently:\n{0}\033[0m".format(vecB)
    print( "\033[0m" , end='' )
    # compute sum
    scalarC = 0.0
    for j in range(3): scalarC += vecA[j] * vecB[j]
    # return
    return scalarC

def vecCross( vecA , vecB ):
    # check the inputs
    print( "\033[3;38;5;235m" , end='' )
    assert isinstance(vecA,list) and len(vecA)==3 and type(vecA[0]) in [int,float,complex], "\n\033[0;38;5;202mERROR: cross() input 1 must be a length-3 list of numbers, currently:\n{0}\033[0m".format(vecA)
    assert isinstance(vecB,list) and len(vecB)==3 and type(vecB[0]) in [int,float,complex], "\n\033[0;38;5;202mERROR: cross() input 2 must be a length-3 list of numbers, currently:\n{0}\033[0m".format(vecB)
    print( "\033[0m" , end='' )
    # initialize result
    vecC = [ 0.0 for i in range(3) ]
    # compute elements
    vecC[0] = vecA[1]*vecB[2] - vecA[2]*vecB[1]
    vecC[1] = vecA[2]*vecB[0] - vecA[0]*vecB[2]
    vecC[2] = vecA[0]*vecB[1] - vecA[1]*vecB[0]
    # return
    return vecC

def vecNorm( vecA ):
    # return norm via dot product
    return vecDot(vecA,vecA)**0.5

def vecReciprocals( vecA1 , vecA2 , vecA3 ):
    # check the inputs
    print( "\033[3;38;5;235m" , end='' )
    assert isinstance(vecA1,list) and len(vecA1)==3 and type(vecA1[0]) in [int,float,complex], "\n\033[0;38;5;202mERROR: cross() input A1 must be a length-3 list of numbers, currently:\n{0}\033[0m".format(vecA1)
    assert isinstance(vecA2,list) and len(vecA1)==3 and type(vecA1[0]) in [int,float,complex], "\n\033[0;38;5;202mERROR: cross() input A1 must be a length-3 list of numbers, currently:\n{0}\033[0m".format(vecA2)
    assert isinstance(vecA3,list) and len(vecA1)==3 and type(vecA1[0]) in [int,float,complex], "\n\033[0;38;5;202mERROR: cross() input A1 must be a length-3 list of numbers, currently:\n{0}\033[0m".format(vecA3)
    print( "\033[0m" , end='' )
    # compute vectors
    (sA1,sA2,sA3) = ( vecNorm(vecA1) , vecNorm(vecA2) , vecNorm(vecA3) )
    (vecB1,vecB2,vecB3) = ( vecCross(vecA2,vecA3) , vecCross(vecA3,vecA1) , vecCross(vecA1,vecA2) )
    (sB1,sB2,sB3) = ( vecNorm(vecB1) , vecNorm(vecB2) , vecNorm(vecB3) )
    # compute normalizing vactors
    sAB1 = vecDot(vecA1,vecB1)
    sAB2 = vecDot(vecA2,vecB2)
    sAB3 = vecDot(vecA3,vecB3)
    # normalize vectors ( up to 2pi )
    vecB1 = [ 2*3.141592653589793/sAB1 * term for term in vecB1 ]
    vecB2 = [ 2*3.141592653589793/sAB2 * term for term in vecB2 ]
    vecB3 = [ 2*3.141592653589793/sAB3 * term for term in vecB3 ]
    # return
    return (vecB1,vecB2,vecB3)





def isKpointListing(inList):
    assert isinstance(inList,list), "isKpointListing() input must be a list!"
    assert len(inList)==3, "isKpointListing() input must have length 3!"
    if inList[2]==0.0 and not any([ abs(val)>0.5 or abs(val)<1E-2 for val in inList[:-1] ]): return True
    return False







# ---------------------------------------------------------------------------- #
if __name__ == '__main__':

    # import commandline arguments
    from vtySH_arguments import arguments
    myArguments = arguments()

    # initialize WAVECAR dictionary
    wavecar = { key:None for key in _keyListWavecar() }

    # append to wavecar
    if not myArguments["-k"]==None: wavecar["k"] = myArguments["-k"]
    if not myArguments["-b"]==None: wavecar["b"] = myArguments["-b"]

    loadWAVECAR("WAVECAR-test-bilayer-ncl/WAVECAR")





# -----------------------------------------------------------------------------#
# -----------------------------------------------------------------------------#
# -----------------------------------------------------------------------------#
__module__    = "VTY"
__status__    = "Prototype" # (Prototype,Development,Production)
__author__    = "(Jonathan) Tyler Reichanadter"
__email__     = "jtreichanadter@berkeley.edu"
__copyright__ = "Copyright 2020, Neaton Group at UC Berkeley"
__license__   = "GPL"
