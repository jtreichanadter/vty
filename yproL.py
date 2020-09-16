#! /usr/bin/env python





# --- IMPORTS ------------------------------------------------------------------



# --- parsing
from _parsing.parsePOSCAR import parsePoscar
from _parsing.parsePROCAR import parseProcar
from _parsing.parseOUTCAR import parseOutcar

# --- bash
from _utilities.shellGuides import BASH

# --- math
from math import pi
from math import atan2 as arctan

# --- variable access
from _utilities.shellGuides import pullArguments
from _utilities.shellGuides import identifyFileName
vars = pullArguments()
fileName = identifyFileName(vars)
from _utilities.shellGuides import string2List
from _utilities.shellGuides import castListInt
from _utilities.shellGuides import castListUint

# --- create escape color strings
from _utilities.shellGuides import bashColor as Color
noColor = Color("reset")
Red     = Color("red",   0)
Blue    = Color("cyan",  0)
Gray    = Color("gray",  1)
Dark    = Color("black", 1)
Green   = Color("green", 1)
Pink    = Color("pink",  1)
Orange  = Color("orange",1)
Purple  = Color("purple",1)
Colors = [ Green , Pink , Orange , Purple , Gray , Red , Blue ]





# --- MAIN ---------------------------------------------------------------------


# --- parse data
try: PROCAR = parseProcar( fileName )
except: PROCAR = parseProcar( identifyFileName(["."]) )
try: OUTCAR = parseOutcar( fileName )
except: OUTCAR = parseOutcar( identifyFileName(["."]) )
try: POSCAR = parsePoscar( fileName )
except: POSCAR = parsePoscar( identifyFileName(["."]) )


# --- initialize arguments
# k-point
K = PROCAR["k"]
k = [1]
for j in range(2,len(K)-1):
    if ( (K[j][0]-K[j-1][0])**2 + (K[j][1]-K[j-1][1])**2 + (K[j][2]-K[j-1][2])**2 ) < 1e-12: k.append(j)
k.append(len(K))
# b-band
B = PROCAR["bands"]
b = []
filled = False
for j in range(len(B[0])):
    if   not filled and PROCAR["occupancies"][0][j]<0.9 and b==[]:
            b = [j-1,j,j+1]
    elif not filled and PROCAR["occupancies"][0][j]>0.1 and not b==[]:
        b.append(j+1)
    elif not filled and not b==[]:
        b.append(j+1)
        b.append(j+2)
        b.append(j+3)
        filled = True
# --- ion list
i = range(len(POSCAR["elements"]))
I = POSCAR["elements"]
# make color array
ColorIons = []
Elems = []
for ii in i:
    if I[ii] not in Elems:
        Elems.append(I[ii])
    ColorIons.append( Colors[ Elems.index(I[ii])] )

# --- interpert variable input
spdBasis = True
if not len(vars)==0:
    for j in range(len(vars)-1):
        if   vars[j] in ["ylm","Ylm","yLM","YLM"]: spdBasis = False
        elif vars[j] in ["-k","-K"]: k = castListUint(string2List( vars[j+1] ))
        elif vars[j] in ["-b","-B"]: b = castListUint(string2List( vars[j+1] ))
# check range of values
assert max(k)<=len(K),    Red + "Error! Largest Kpt index is {0}".format(len(K)) + noColor
assert min(k)>0,          Red + "Error! Kpt indices must be positive!" + noColor
assert max(b)<=len(B[0]), Red + "Error! Largest band index is {0}".format(len(B[0])) + noColor
assert min(b)>0,          Red + "Error! Band indices must be positive!" + noColor


# --- acquire PROCAR data
labels = PROCAR["projectionLabels"]
proj  = PROCAR["projections"]
projC = PROCAR["projectionsC"]
Lq = range(len(PROCAR["projections"]))
Ln = range(len(PROCAR["projections"][0]))
Li = range(len(PROCAR["projections"][0][0]))
Ls = range(len(PROCAR["projections"][0][0][0]))
projT = [[ PROCAR["projectionsTotal"][q][n][-1] for n in Ln] for q in Lq]
muZ = PROCAR["muProjectionZ"]

# --- create printing header
headerString = "KPT BAND | (E-Ef) (occ) |   Psi "
if spdBasis: headerString += "(spd)"
else:        headerString += "(ylm)"

# --- acquiring phi of orbital
def getAngleStr( ValC , *offset ):
    # grab optional argument
    try: offset = offset[0]
    except: offset = 0.0
    # initalize angle array
    angle = arctan( ValC.imag , ValC.real ) - offset
    #print("Angle = {0}".format(angle))
    if   angle < -pi: angle = angle + pi
    elif angle >  pi: angle = angle - pi
    tempList = [ abs(angle-j*pi/6) for j in range(-6,7) ]
    idx = tempList.index( min(tempList) )
    if    idx == 0: # phi ~ -pi
        return " +" + Dark + "e^(-6/6pi)  " + noColor
        return " - "
    elif idx == 1: # phi ~ (-5/6)pi
        return " +" + Dark + "e^(-5/6pi)  " + noColor
        return " -" + Dark + "e^(pi/6)  " + noColor
    elif idx == 2: # phi ~ (-4/6)pi
        return " +" + Dark + "e^(-4/6pi)  " + noColor
        return " -" + Dark + "e^(pi/3)  " + noColor
    elif idx == 3: # phi ~ (-3/6)pi
        return " +" + Dark + "e^(-3/6pi)  " + noColor
        return " -" + Dark + "i " + noColor
    elif idx == 4: # phi ~ (-2/6)pi
        return " +" + Dark + "e^(-2/6pi)  " + noColor
        return " -" + Dark + "e^(2pi/3) " + noColor
    elif idx == 5: # phi ~ (-1/6)pi
        return " +" + Dark + "e^(-1/6pi)  " + noColor
        return " -" + Dark + "e^(5pi/6) " + noColor
    elif idx == 6: # phi ~ 0
        return " +" + Dark + "e^(00/6pi)  " + noColor
        return " + "
    elif idx == 7: # phi ~ (+1/6)pi
        return " +" + Dark + "e^(+1/6pi)  " + noColor
        return " +" + Dark + "e^(pi/6)  " + noColor
    elif idx == 8: # phi ~ (+2/6)pi
        return " +" + Dark + "e^(+2/6pi)  " + noColor
        return " +" + Dark + "e^(pi/3)  " + noColor
    elif idx == 9: # phi ~ (+3/6)pi
        return " +" + Dark + "e^(+3/6pi)  " + noColor
        return " +" + Dark +  "i " + noColor
    elif idx == 10: # phi ~ (+4/6)pi
        return " +" + Dark + "e^(+4/6pi)  " + noColor
        return " +" + Dark + "e^(2pi/3) " + noColor
    elif idx == 11: # phi ~ (+5/6)pi
        return " +" + Dark + "e^(+5/6pi)  " + noColor
        return " +" + Dark + "e^(5pi/6) " + noColor
    elif idx == 12: # phi ~ (+6/6)pi
        return " +" + Dark + "e^(+6/6pi)  " + noColor
        return " - "

# --- print
# print k-vector header
for kk in k:
    # form dark vector
    v1 = str(round(K[kk-1][0],2)).rjust(5)
    v2 = str(round(K[kk-1][1],2)).rjust(5)
    v3 = str(round(K[kk-1][2],2)).rjust(5)

# print elements header
print("--> Elements")
for ii in i:
    posVector = " "
    for j in range(3): posVector = posVector + str(round(POSCAR["positionsDirect"][ii][j],2)).ljust(4,"0") + " "
    print("{0}{1} = {2} ({3}){4}".format( ColorIons[ii] , str(ii+1).rjust(2) , (str(I[ii]).strip()).ljust(2) , posVector , noColor ))
# print LOOP
print("--> Projections")
for kk in k:
    print(headerString)
    for bb in b:
        # --- CONSTRUCT WAVE-FUNCTION
        WAVE = Dark + "(" + str(round(projT[kk-1][bb-1]*100)) + "%) " + noColor
        for i in Li:
            for s in range(len(Ls)-1):
                # check for significant contribution:
                if proj[kk-1][bb-1][i][s] >= 0.1:
                    # L-phi
                    #try:
                    #    phiStr = getAngleStr( projC[kk-1][bb-1][i][s] , offset )
                    #except:
                    #    phiStr = " "
                    #    offset = arctan( projC[kk-1][bb-1][i][s].imag , projC[kk-1][bb-1][i][s].real )
                    phiStr = getAngleStr( projC[kk-1][bb-1][i][s] )
                    # mu-Z
                    muValZ = muZ[kk-1][bb-1][i][s] / proj[kk-1][bb-1][i][s]
                    if muValZ<=-0.1: muZstr = Blue + "(" + str(round(muValZ,2)).ljust(5,"0")  + ")" + noColor
                    elif muValZ<0.0: muZstr = Dark + "(" + str(round(muValZ,2)).ljust(5,"0")  + ")" + noColor
                    elif muValZ<0.1: muZstr = Dark + "(+" + str(round(muValZ,2)).ljust(4,"0") + ")" + noColor
                    else:            muZstr = Red  + "(+" + str(round(muValZ,2)).ljust(4,"0") + ")" + noColor
                    # append to wave string
                    WAVE = WAVE + phiStr + str(round(proj[kk-1][bb-1][i][s],2)).ljust(4,"0") + " " + ColorIons[i] + str(i+1) + labels[s] + muZstr + " "
        WAVE = WAVE.strip("+").strip()
        # finalize variables
        kIdx = str(kk).rjust(3)
        bIdx = str(bb).rjust(3)
        Ef   = round(OUTCAR["bandsFermi"][kk-1][bb-1],3)
        Occ  = str(round(PROCAR["occupancies"][kk-1][bb-1],2)).ljust(4)
        if Ef>0: Ef  = " " + str(Ef).ljust(6)
        else:    Ef = str(Ef).ljust(7)
        if PROCAR["occupancies"][kk-1][bb-1]<0.5: Occ = Dark + Occ + noColor
        # print lines
        print("{0} {1}  | {2} {3} |  {4}".format(kIdx,bIdx,Ef,Occ,WAVE)  )
        #del offset










#
