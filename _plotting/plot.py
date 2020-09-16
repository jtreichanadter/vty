#! /usr/bin/env python

# --- IMPORTS ------------------------------------------------------------------

# ### bash ###
from subprocess import call as spCall
def BASH(command):
    if isinstance(command,list): spCall(command)
    elif isinstance(command,str): spCall(command.split())
    else: warning("could not interperet bash command")
#def BASHQ(command):
#    try:
#        if isinstance(command,str): command = command.split()
#        if isinstance(command,list): spCall( command + [">",".tmp.trashi"] )
#        else: warning("could not interperet bash command")
#    except: pass
# BASH("if test -d '__pycache__'; then rm -r __pycache__; fi")










# ############################# CORE FUNCTIONS ############################### #



# ##### MAIN PLOTTING UTILITY ##### #
def plot(x,y,*properties):
# input:
#  x - Xaxis data, 1D or 1Dx3 array
#       if 1D, assume simple plot
#       if 2D, assume k-vector input and construct paths
#  y - Yaxis data, nD array, 1st dimension must match 'x'
#       if 1D, assume simple plot
#       if 2D, assume bandstructure, although the
# properties - optional dictionary to specify plot details

    # ----- verify data input
    assert isinstance(x,list), error( "input 'x' to plot() was not a list :(","x = {0}".format(x) )
    assert isinstance(y,list), error( "input 'y' to plot() was not a list :(","y = {0}".format(y) )
    if isinstance(y[0],list): assert not isinstance(y[0][0],list), error( "input 'y' to plot() exceeds 2 dimensions!","size(x) = {0} x {1} x {2} x ...".format(L1,L2,len(y[0][0])) )
    assert len(x)==len(y), error( "inputs 'x','y' to plot() are not equal length :(","len(x)={0}  len(y)={1}".format(len(x),len(y)) )
    L1 = len(x)
    if isinstance(y[0],list): L2 = len(y[0])
    else: L2 =1

    # ----- manage 'properties' input
    try: properties = properties[0]
    except: properties = makeDictionaryProperties()

    # ----- format data and 'properties' for gnu-plot
    ( x , y , properties ) = groomProperties( x , y , properties )
    fileName = properties["fileName"]
    xMin = properties["rangeX"][0]
    xMax = properties["rangeX"][1]
    yMin = properties["rangeY"][0]
    yMax = properties["rangeY"][1]
    borderThickness = properties["borderThickness"]
    thickness = properties["thickness"]
    titleFontSize = properties["titleFontSize"]
    title = properties["title"]
    labelX = properties["labelX"]
    labelY = properties["labelY"]
    ticsX = properties["ticsX"]
    color = properties["color"]
    aspectRatio = properties["aspectRatio"]

    # ----- write data file
    if isinstance(y[0],list): open(".temporary.dat","w").write("\n".join("{0} {1}".format(x[i]," ".join(str(y[i][j]) for j in range(L2)  ) ) for i in range(L1) ))
    else: open(".temporary.dat","w").write("\n".join("{0} {1}".format(x[i],y[i]) for i in range(L1) ))

    # ----- write meta for gnu plot script
    GNU  = open(".temporary.gnu","w")
    GNU.write("# ########################################### #\n")
    GNU.write("# ##### PYTHON GENERATED GNUPLOT SCRIPT ##### #\n")
    GNU.write("# ########################################### #\n")
    # [0] set output figure details (color '.eps' , fileName)
    GNU.write("\n# --- 0 --- set output figure format to a color '.eps' file\n")
    GNU.write( "set terminal postscript color\n" )
    GNU.write( "set termoption enhanced\n" )
    GNU.write( "set encoding utf8\n" )
    GNU.write( "set output '{0}'\n".format(fileName) )
    GNU.write( "set size ratio {0}\n".format(aspectRatio) )
    GNU.write( "set multiplot\n" )
    # [1] border and axes details
    GNU.write( "\n# --- 1 --- establish plotting ranges and size\n" )
    # (1a) set x,y ranges
    GNU.write( "set xrange [ {0} : {1} ]\n".format(xMin,xMax) )
    GNU.write( "set yrange [ {0} : {1} ]\n".format(yMin,yMax) )
    GNU.write( "set border lw {0}\n\n".format(borderThickness) )
    # [2] specify legends, keys, and add grid
    GNU.write( "\n# --- 2 --- specify legends, keys, and include grid\n" )
    GNU.write( "unset colorbox\n" )
    GNU.write( "set key off\n" )
    GNU.write( "set grid\n" )
    # [3] set plot Labels
    GNU.write( "\n# --- 3 --- set plot labels\n" )
    # (3a) labels
    GNU.write( "set title '{0}'\n".format(title) )
    GNU.write( "set xlabel '{0}'\n".format(labelX) )
    GNU.write( "set ylabel '{0}'\n".format(labelY) )
    # (3v) font sizes
    if titleFontSize!=None:
        GNU.write(  "set title font ',{0}'\n".format(titleFontSize   ) )
        GNU.write( "set ylabel font ',{0}'\n".format(titleFontSize-6 ) )
        GNU.write(  "set xtics font ',{0}'\n".format(titleFontSize-6 ) )
        GNU.write(  "set ytics font ',{0}'\n".format(titleFontSize-10) )
    # [4] set x-tics
    GNU.write( "set xtics {0}\n".format(ticsX) )
    #GNU.write( "set xtics {0},1,{1}\n".format(min(x),max(x)) )
    #for i in range(len(ticsX)):
        #GNU.write( "set xtics add ('{1}' {0})\n".format(ticsX[i][0],ticsX[i][1]) )
    GNU.write( "set xtics scale 0,1\n" )


    # ----- write operational code for gnu plot script
    GNU.write( "\n# --- 5 --- specify plot styles\n" )
    if not isinstance(y[0],list):
        # -- CASE 1: Y=1D , color=0D
        if not isinstance(color,list):
            # set style
            GNU.write( "set style line 1 lc rgb '{0}' lw {1}\n".format(color,thickness[0]) )
            # plot
            GNU.write( "\n# --- 6 --- PLOT\n" )
            if properties["draw0X"]: GNU.write( "plot 0 w lines lc rgb '#000000' lw 3\n" )
            GNU.write( "plot '.temporary.dat' w lines ls (1)\n" )
        # -- CASE 2: Y=1D , color=1D
        else:
            # set style
            print(len(color))
            for q in range(1,L1):
                GNU.write( "set style line {0} lc rgb '{1}' lw {2}\n".format(q,color[q-1],thickness[0]) )
            # plot
            GNU.write( "\n# --- 6 --- PLOT\n" )
            if properties["draw0X"]: GNU.write( "plot 0 w lines lc rgb '#000000' lw 3\n" )
            GNU.write( "plot for [s=1:{0}] '.temporary.dat' every ::(s)::(s+1) w lines ls (s)\n".format(L1-1) )
    else:
        # -- CASE 3: Y=2D , color=0D
        if not isinstance(color,list):
            # set style
            GNU.write( "set style line 1 lc rgb '{0}' lw {1}\n".format(color,thickness[0]) )
            # plot
            GNU.write( "\n# --- 6 --- PLOT\n" )
            if properties["draw0X"]: GNU.write( "plot 0 w lines lc rgb '#000000' lw 3\n" )
            GNU.write( "plot for [s = 2:{0}] '.temporary.dat' using 1:(column(s)) w lines ls (1)\n".format(L2+1) )
        # -- CASE 4: Y=2D , color=1D
        elif not isinstance(color[0],list):
            # set style
            for b in range(L2):
                GNU.write( "set style line {0} lc rgb '{1}' lw {2}\n".format(b+1,color[b],thickness[0]) )
            # plot
            GNU.write( "\n# --- 6 --- PLOT\n" )
            if properties["draw0X"]: GNU.write( "plot 0 w lines lc rgb '#000000' lw 3\n" )
            GNU.write( "plot for [s = 2:{0}] '.temporary.dat' using 1:(column(s)) w lines ls (s-1)\n".format(L2+1) )
        # -- CASE 5: Y=2D , color=2D
        else:
            # assemble coloring sub-array
            L = (L1-1)*L2
            GNU.write( "\n# -- 5.5 -- plot array\n" )
            GNU.write( "array B[{0}]\n".format(L) )
            GNU.write( "array Q[{0}]\n".format(L) )
            GNU.write( "do for [q=1:{0}] {{\n".format(L1-1) )
            GNU.write( "  do for [b=0:{0}] {{\n".format(L2-1) )
            GNU.write( "     B[q+{0}*b] = b+2\n".format(L1-1) )
            GNU.write( "     Q[q+{0}*b] = q+1\n".format(L1-1) )
            GNU.write( "  }\n" )
            GNU.write( "}\n" )
            # set style
            for q in range(1,L1):
                for b in range(L2):
                    GNU.write( "set style line {0} lc rgb '{1}' lw {2}\n".format(q+(L1-1)*b,color[q-1][b],thickness[q-1][b]) )
            # plot
            GNU.write( "\n# --- 6 --- PLOT\n" )
            if properties["draw0X"]: GNU.write( "plot 0 w lines lc rgb '#000000' lw 3\n" )
            GNU.write( "plot for [s=1:{0}] '.temporary.dat' using 1:(column(B[s])) every ::Q[s]-1::Q[s] w lines ls (s)\n".format(L) )

    # EXECUTE
    GNU.write("\n\n\n")
    GNU.close()
    BASH( "gnuplot .temporary.gnu" )
    BASH( "open {0}".format(fileName) )













# ############################# HELPER FUNCTIONS ############################# #



# prints 'plot' warnings
def warning(*argv):
    print("\033[1;33m--- WARNING : plot.py ---\033[0;33m")
    for arg in argv:
        print(arg)
    print("\033[0m")


# prints 'plot' errors
def error(*argv):
    print("\033[1;31m--- ERROR: plot.py ---\033[0;31m")
    for arg in argv:
        print(arg)
    print("\033[0m")


# avoids file overwriting
def dodgeOverwrite(name):
# checks if 'name' already exists,
# and increments the name until no such file exists.
    # verify string input
    assert isinstance(name,str), error("input 'name={0}' to dodgeOverwrite() was not a string :(".format(name))
    # check base case, return 'name' if it is not a file already
    try:
        tmp = open(name,'r')
        tmp.close()
    except:
        return name
    # modify for file extension
    suffix = ""
    if name.rfind('.')>0:
        suffix = name[name.rfind('.'):]
        name = name[:name.rfind('.')]
    # iterate until a new fileName is formed
    for j in range(1000):
        try:
            NAME = name + "-" + str(j).zfill(3) + suffix
            tmp = open(NAME,'r')
            tmp.close()
            j = j+1
        except:
            return NAME
    # failed to create unique file-name
    warning( "function dodgeOverwrite() FAILED to make a unique file name for '{0}'".format(name) , "new file name is '{0}'".format(name+"SORRY"+suffix) )
    return name+"SORRY"+suffix





# ############################ PROPERTY HANDLING ############################# #



# create properties dictionary
def keyListProperties():
    return [ "aspectRatio" , \
             "borderThickness" , \
             "color" , "colorFormat" , \
             "draw0X" , \
             "fileName" , \
             "labelX" , "labelY" , \
             "rangeX" , "rangeY" , \
             "ticsX" , "ticsY" , \
             "title" , "titleFontSize" , \
             "thickness" ]
# STUB

# verify properties dictionary
def isDictionaryProperties( dictionary ):
    # verifiy dictionary data-type input
    if not isinstance(dictionary,dict):
        error( "Input to isDictionaryProperties() is not even a dictionary data type!","input_type = {0}".format(type(dictionary)),"input = \n{0}".format(dictionary) )
        return False
    # verify input contains every property key
    if any( not key in dictionary.keys() for key in keyListProperties() ):
        tmpSTR = []
        for key in keyListProperties():
            if not key in dictionary.keys(): tmpSTR.append(key)
        error( "Input to isDictionaryProperties() lacks necessary keys.","Missing the following: \n{0}".format(tmpSTR) )
        return False
    # check if input contains any exrta fields
    if any( not key in keyListProperties() for key in dictionary.keys() ):
        tmpSTR = []
        for key in dictionary.keys():
            if not key in keyListProperties(): tmpSTR.append(key)
        warning( "Input to isDictionaryProperties() contains extra keys not part of a properties dictionary.","Extra keys include: \n{0}".format(tmpSTR) )
        return True
    # otherwise return true
    return True


# initialize a properties dictionary
def makeDictionaryProperties():
    return { key:None for key in keyListProperties() }


# print user help guide
def explainProperties( *key ):
# optional input: key - will just print documentation for that particular property key.

    # ### handle 'key' input
    try:
        key = str(key[0])
        if not key in keyListProperties:
            warning("input 'key={0}' to explainProperties() was not recognized :(".format(key))
    except: key = ""

    #
    print("--- PLOT PROPERTIES (python dictionary for gnuplot) ---")
    print(" ")
    print("Tyler's python plotting scripts exploit the gnuplot utility.")
    print("The workflow of this python script is: writeDat(.dat)-->writeGnuScript(.gnu)-->executeGnu(bash)")
    print("To specify plot properties, most python functions accept an optional 'properties' dictionary")
    print("This hopefully explains how to utilize 'properties' in order to make some slick plots!")
    print(" ")
    print(" ")
    print("labelX - string argument, label on the horizontal axis")
    print(" ")
    print("labelXsize - ")
# STUB


# convert properties input for gnu format
def groomProperties( x , y , properties ):

    # [0] check valid dictionary
    assert isDictionaryProperties(properties), error( "input to groomProperties() is not a valid properties dictionary!" )

    # [1] acquire data dimensions
    L1 = len(x)
    if isinstance(y[0],list): L2 = len(y[0])
    else: L2 = 1


    # [2] modify properties into gnuplot format

    # ---- rangeX
    X = [ x[i] for i in range(L1) ]
    key = "rangeX"
    value = properties[key]
    if value!=None:
        # validate input
        assert isinstance(value,list), error( "'{0}' property is not a list!".format(key),"properties['{0}'] = \n{1}".format(key,value) )
        assert len(value)==2, error( "'{0}' property is not length 2".format(key),"properties['{0}'] = \n{1}".format(key,value) )
        assert isinstance(value[0],(int,float)), error( "'{0}' property has inappropriate entries".format(key),"properties['{0}'].type =\n{1}".format(key,type(value[0])) )
        assert value[0]<value[1], error( "'{0}' property is not a logically ordered range".format(key),"properties['{0}'].type ={1}".format(key,value) )
        assert value[0]<value[1], error( "'{0}' property is not a logically ordered range".format(key),"properties['{0}'].type ={1}".format(key,value) )
    else:
        if not isinstance(x[0],list):
            # assign values based on 'x' data
            properties[key] = [0,0]
            properties[key][0] = X[0]  - 0.001 * ( (X[-1]-X[0]) + 1e-8 )
            properties[key][1] = X[-1] + 0.001 * ( (X[-1]-X[0]) + 1e-8 )
            xMin = properties[key][0]
            xMax = properties[key][1]
        else:
            properties[key] = [0,0]
            properties[key][0] = 0
            properties[key][1] = 1
            xMin = 0
            xMax = 1

    # ---- modify 'x' data (and ticsX) as needed for k-vectors
    # handle k-vector input
    if isinstance(x[0],list):
        # verify 3-vector input
        assert len(x[0])==3, error( "input 'x' to plot() was multidimensional,","but the 2nd dimension length is not 3 (for k-vector input)","size(x) = {0} x {1}".format(L,len(x[0])) )
        assert not isinstance(x[0][0],list), error( "input 'x' to plot() exceeds 2 dimensions!","size(x) = {0} x {1} x {2} x ...".format(L,3,len(x[0][0])) )
        # convert to 1D array
        X[0] = 0.0
        dk = [ ( (x[i][0]-x[i-1][0])**2 + (x[i][1]-x[i-1][1])**2 + (x[i][2]-x[i-1][2])**2 )**(0.5) for i in range(1,L1) ]
        dk_avg = sum(dk)/(L1-1)
        for i in range(1,L1):
            if dk[i-1] < 0.001*dk_avg: dk[i-1] = 0
            X[i] = X[i-1] + dk[i-1]
        # normalize 1D array
        for i in range(1,L1): X[i] = X[i]/X[-1]
        # ---- rangeX (modify for k-path input)
        if properties["rangeX"][1]>1:
            properties["rangeX"][0] = X[ properties["rangeX"][0] ] - 1e-8
            properties["rangeX"][1] = X[ properties["rangeX"][1] ] + 1e-8
            xMin = properties["rangeX"][0]
            xMax = properties["rangeX"][1]
        # ---- ticksX
        ticsX = properties["ticsX"]
        #   case 1: proper input, list of len2 lists [ point , label ]
        if isinstance(ticsX,list) and isinstance(ticsX[0],list) and len(ticsX[0])==2:
            kIdxFormat = False
            for i in range( len(ticsX)-1 , -1 , -1 ):
                if kIdxFormat or ticsX[i][0]>1:
                    ticsX[i][0] = X[ int(round(ticsX[i][0])) ]
                    kIdxFormat = True
        #   case 2: no input, identify ticks by finding when dk=0
        elif ticsX==None:
            ticsX = []
            for i in range(1,L1):
                tmpSTR = "(" + str(round(x[i][0],2)) + "," + str(round(x[i][1],2)) + "," + str(round(x[i][2],2)) + ")"
                if dk[i-1]==0: ticsX.append( [ X[i] , tmpSTR ] )
        #   case 3: random input, identify ticks by finding when dk=0, no label
        else:
            ticsX = []
            for i in range(1,L1):
                if dk[i-1]==0: ticsX.append( [ X[i] , " " ] )
    # ordinary case ticks handling
    else:
        # ---- ticsX
        ticsX = properties["ticsX"]
        key = "ticsX"
        value = properties[key]
        if value!=None:
            # case 1: a perfect nested list
            if isinstance(ticsX,list) and isinstance(ticsX[0],list) and len(ticsX[0])==2:
                for i in range(len(value)):
                    assert isinstance(ticsX[i][0],(float,int)), error( "'{0}' property is not a nested list with [ [ point(float) , label(str)  ] , ... ] form!".format(key),"properties['{0}'] = \n{1}".format(key,value) )
                    assert isinstance(ticsX[i][1],str), error( "'{0}' property is not a nested list with [ [ point(float) , label(str)  ] , ... ] form!".format(key),"properties['{0}'] = \n{1}".format(key,value) )
            # case 2: list of points
            elif isinstance(ticsX,list) and isinstance(ticsX[0],(int,float)):
                for i in range(len(value)-1,-1,-1):
                    if ticsX[i]<properties["rangeX"][0] or ticsX[i]>properties["rangeX"][1]:
                        del ticsX[i]
                    else:
                        ticsX[i] = [ ticsX[i] , str(ticsX[i]) ]
            # case 3: single value for subdivision
            assert isinstance(ticsX,(int,float)), error( "'{0}' property is not a recognized form!".format(key),"properties['{0}'] = \n{1}".format(key,value) )
            Lt = ticsX+2
            ticsX = [ [0,"0"] for i in range(Lt) ]
            for i in range(Lt):
                j =  xMin + ( i / (Lt-1) ) * (xMax-xMin)
                ticsX[i] = [ j , str(round(j+1e-8,1)) ]
    # FINALLY: formalize x-ticks into a single string
    tmpSTR = "("
    for j in range(len(ticsX)):
        tmpSTR = tmpSTR + "'{1}' {0}, ".format(ticsX[j][0],ticsX[j][1])
    tmpSTR = tmpSTR.rstrip(", ") + ")"
    properties["ticsX"] = tmpSTR

    # ---- rangeY
    key = "rangeY"
    value = properties[key]
    if value!=None:
        # validate input
        assert isinstance(value,list), error( "'{0}' property is not a list!".format(key),"properties['{0}'] = \n{1}".format(key,value) )
        assert len(value)==2, error( "'{0}' property is not length 2".format(key),"properties['{0}'] = \n{1}".format(key,value) )
        assert isinstance(value[0],(int,float)), error( "'{0}' property has inappropriate entries".format(key),"properties['{0}'].type =\n{1}".format(key,type(value[0])) )
    else:
        if isinstance(y[0],list):
            properties[key] = [0,0]
            properties[key][0] = min( [min(y[i]) for i in range(L1)] )
            properties[key][1] = max( [max(y[i]) for i in range(L1)] )
        else:
            yMin = min( y )
            yMax = max( y )
            properties[key] = [0,0]
            properties[key][0] = yMin - 0.1 * ( (yMax-yMin) + 1e-8 )
            properties[key][1] = yMax + 0.1 * ( (yMax-yMin) + 1e-8 )

    # ---- borderThickness
    key = "borderThickness"
    value = properties[key]
    if not properties[key]==None:
        assert isinstance(value,(int,float)), error( "'{0}' property has inappropriate entries".format(key),"properties['{0}'].type =\n{1}".format(key,type(value)) )
    else:
        properties[key] = 1.5

    # ---- aspectRatio
    key = "aspectRatio"
    value = properties[key]
    if not properties[key]==None:
        assert isinstance(value,(float,int)), error( "'{0}' property has an inappropriate entry".format(key),"properties['{0}'].type =\n{1}".format(key,type(value)) )
    else:
        properties[key] = 0.5

    # ---- labelX
    key = "labelX"
    value = properties[key]
    if not properties[key]==None:
        assert isinstance(value,str), error( "'{0}' property has inappropriate entries".format(key),"properties['{0}'].type =\n{1}".format(key,type(value)) )
    else:
        properties[key] = "X"

    # ---- labelY
    key = "labelY"
    value = properties[key]
    if not properties[key]==None:
        assert isinstance(value,str), error( "'{0}' property has inappropriate entries".format(key),"properties['{0}'].type =\n{1}".format(key,type(value)) )
    else:
        properties[key] = "Y"

    # ---- titleFontSize
    key = "titleFontSize"
    value = properties[key]
    if not properties[key]==None:
        assert isinstance(value,(int,float)), error( "'{0}' property has inappropriate entries".format(key),"properties['{0}'].type =\n{1}".format(key,type(value)) )

    # ---- title
    key = "title"
    value = properties[key]
    if not properties[key]==None:
        assert isinstance(value,str), error( "'{0}' property has inappropriate entries".format(key),"properties['{0}'].type =\n{1}".format(key,type(value)) )
    else:
        properties[key] = " "

    # ---- fileName
    key = "fileName"
    value = properties[key]
    if properties[key]==None:   properties[key] = dodgeOverwrite(".temporary.eps")
    else:
        if properties[key][-4:]==".eps": fileName = properties[key]
        else: fileName = properties[key]+".eps"
        if fileName!=dodgeOverwrite(fileName): warning( "function dodgeOverwrite() found '{0}' was taken.".format(fileName) , "Instead, file will be named '{0}'".format(dodgeOverwrite(fileName)) )
        properties[key] = dodgeOverwrite(fileName)

    # ---- draw0X
    key = "draw0X"
    value = properties["draw0X"]
    if not properties[key]==None and not properties[key]:
        properties[key] = False
    else:
        properties[key] = True

    # ---- color
    key = "color"
    value = properties[key]
    if value!=None:
        # access color input type
        colorCase = checkColorType( value )
        # case 01: single string, name
        if colorCase==1:
            properties[key] = colorName2Hex(value)
        # case 02: single 3vector (rgb,hsv)
        elif colorCase==2:
            properties[key] = colorRGB2Hex(value)
        # case 11: 1D list of strings, names
        elif colorCase==11:
            if L2==1: L = L1
            else: L = L2
            assert L==len(value), error("color array for plot() is insuffient size!.","There are dim(y)={0}x{1} plots and {2} colors".format(L1,L2,len(value)) )
            properties[key] = [ colorName2Hex(value[i]) for i in range(L) ]
        # case 12: 1D list, list of 3vectors (rgb,hsv)
        elif colorCase==12:
            if L2==1:
                L = L1
                assert L==len(value), error("color array for plot() is insuffient size!.","There are dim(y)={0}x{1} plots and {2} colors".format(L1,L2,len(value)) )
                # make simple interpolate
                properties[key] = [ colorRGB2Hex([ (value[i][0]+value[i-1][0])/2 , (value[i][1]+value[i-1][1])/2 , (value[i][2]+value[i-1][2])/2 ]) for i in range(1,len(value)) ]
            else:
                L = L2
                assert L==len(value), error("color array for plot() is insuffient size!.","There are dim(y)={0}x{1} plots and {2} colors".format(L1,L2,len(value)) )
                # convert values
                properties[key] = [ colorRGB2Hex(value[i]) for i in range(L) ]
        # case 21: 2D list of strings, names
        elif colorCase==21:
            assert L1==len(value) and L2==len(value[0]), error("color array for plot() is insuffient size!.","There are dim(y)={0}x{1} plots and colors={2}x{3}".format(L1,L2,len(value),len(value[0])) )
            properties[key] = [ [ colorName2Hex(value[i][j]) for j in range(len(value[0])) ] for i in range(len(value)) ]
        # case 22: 2D list, list of 3vectors (rgb,hsv)
        elif colorCase==22:
            assert L1==len(value) and L2==len(value[0]), error("color array for plot() is insuffient size!.","There are dim(y)={0}x{1} plots and colors={2}x{3}".format(L1,L2,len(value),len(value[0])) )
            properties[key] = [ [ colorRGB2Hex([ (value[i][j][0]+value[i-1][j][0])/2 , (value[i][j][1]+value[i-1][j][1])/2 , (value[i][j][2]+value[i-1][j][2])/2 ]) for j in range(len(value[0])) ] for i in range(1,len(value)) ]
    else:
        properties[key] = "#666666"

    # ---- thickness
    key = "thickness"
    value = properties[key]
    if not properties[key]==None:
        assert isinstance(value,(int,float,list)), error( "'{0}' property has inappropriate entries".format(key),"properties['{0}'].type =\n{1}".format(key,type(value)) )
        tmp = properties["color"]
        if isinstance(value,(int,float)) and isinstance(tmp,list) and isinstance(tmp[0],list):
            properties[key] = [ [ value for j2 in range(len(tmp[0])) ] for j1 in range(len(tmp)) ]
    else:
        tmp = properties["color"]
        if isinstance(tmp,list) and isinstance(tmp[0],list):
            properties[key] = [ [ 1.5 for j2 in range(len(tmp[0])) ] for j1 in range(len(tmp)) ]
        else:
            properties[key] = 1.5




    return ( X , y, properties )
# STUB



# ------------ color-specific functions for type casting --------------


def checkColorType( color ):

    # input is a list ( 0D rgb or 1D+ array )
    if isinstance(color,list):
        # check 0D rgb case
        if len(color)==3 and not isinstance(color[0],(str,list)): return 2
        # input is a nested list ( 1D rgb or 2D array )
        if isinstance(color[0],list):
            # check 1D rgb case
            if len(color[0])==3 and not isinstance(color[0][0],(str,list)): return 12
            # input is a 3D list, must be a 2D rgb
            if isinstance(color[0][0],list):
                assert len(color[0][0])==3 and not isinstance(color[0][0][0],(str,list)), error( "Input 'color' property for plot not recognized, 3D array not rgb.","color = \n{0}".format(color) )
                return 22
            # otherwise must be 2D string case
            assert isinstance(color[0][0],str), error( "Input 'color' property for plot not recognized, 2D array not rgb or string.","color = \n{0}".format(color) )
            return 21
        # otherwise must be 1D string case
        assert isinstance(color[0],str), error( "Input 'color' property for plot not recognized, 1D array not rgb or string.","color = \n{0}".format(color) )
        return 11

    # otherwise must be 0D string case
    assert isinstance(color,str),   error( "Input 'color' property for plot not recognized, not rgb or string.","color = {0}".format(color) )
    return 1


def hsv2rgb():
    pass
# STUB


def colorRGB2Hex( colorRGB ):
    # verify input
    assert isinstance(colorRGB,list), error( "input 'colorRGB' into colorRGB2Hex() is not a list!","colorRGB = {0}".format(colorRGB) )
    assert len(colorRGB)==3, error( "input 'colorRGB' into colorRGB2Hex() is not length-3!","colorRGB = \n{0}".format(colorRGB) )
    # map colors
    for j in range(3):
        colorRGB[j] = abs(colorRGB[j])
        if colorRGB[j]<=1: colorRGB[j] = round(255*colorRGB[j])
        else: colorRGB[j] = min(round(colorRGB[j]),255)
    (R,G,B) = colorRGB
    #cast hex string
    return "#{0}{1}{2}".format( '{:02x}'.format(R) , '{:02x}'.format(G) , '{:02x}'.format(B) )


def colorName2Hex( colorName ):
# ( pre-defined color names follow matlab 2019 convention )
    # verify input
    assert isinstance(colorName,str), error( "input 'colorName' into colorName2Hex() is not a string!","colorName = \n{0}".format(colorName) )
    colorName = colorName.strip()
    # first check if it's already a hex value
    if colorName[0]=="#" and len(colorName)==7: return colorName
    # assemble color dictionary (end of each list must be the rgb-hex value)
    colorDictionary = {}
    colorDictionary["yellow"]  = [ "#ffff00" , "y" , "Y" , "yellow"  , "Yellow"  , "YELLOW"  ]
    colorDictionary["magenta"] = [ "#ff00ff" , "m" , "M" , "magenta" , "Magenta" , "MAGENTA" ]
    colorDictionary["cyan"]    = [ "#00ffff" , "c" , "C" , "cyan"    , "Cyan"    , "CYAN"    ]
    colorDictionary["red"]     = [ "#ff0000" , "r" , "R" , "red"     , "Red"     , "RED"     ]
    colorDictionary["green"]   = [ "#00ff00" , "g" , "G" , "green"   , "Green"   , "GREEN"   ]
    colorDictionary["blue"]    = [ "#0000ff" , "b" , "B" , "blue"    , "Blue"    , "BLUE"    ]
    colorDictionary["white"]   = [ "#ffffff" , "w" , "W" , "white"   , "White"   , "WHITE"   ]
    colorDictionary["black"]   = [ "#000000" , "k" , "K" , "black"   , "Black"   , "BLACK"   ]
    # cross-reference colors with 'colorName'
    for key,colors in colorDictionary.items():
        if colorName in colors: return colors[0]
    # otherwise accept failure
    warning( "Input 'colorName' into colorName2Hex() cannot be identified and mapped to hex" , " colorName = '{0}'\n".format(colorName) , "Here are the acceptable pre-defined colors:\n{0}".format("\n".join("{0} - {1} ".format(key , colorDictionary[key]) for key in colorDictionary )) , "\nDefaulting to color = #666666" )
    return "#666666"






# ################################# RUN CODE ################################# #
if __name__ == '__main__':

    # create temp data
    #x = [ j/10 for j in range(-50,50+1) ]
    x = [ [0,0,0] , [0.2,0.2,0.2] , [0.2,0.2,0.2] , [0.4,0.4,0.4] , [0.5,0.5,0.5] , [0.6,0.6,0.6] , [0.6,0.6,0.6] , [1,1,1] ]
    #y = [ [ i/2+(j/10)**2 for i in range(3)] for j in range(-50,50+1) ]
    #y = [ 1+(j/10)**3 for j in range(-50,50+1) ]
    y = [ 1 , 2 , 2 , 3 , 4 , 5 , 5 , 6]

    # set plot properties
    properties = makeDictionaryProperties()
    #properties["rangeX"] = [0,3]
    #properties["draw0X"] = False
    #properties["color"] = "cool"
    #properties["fileName"] = ".myFigure.eps"
    properties["ticsX"] = 3
    properties["labelX"] = " "
    properties["labelY"] = "energy (E-E_f ,   eV)"
    properties["title"] = "Bandstructure"
    properties["titleFontSize"] = 24

    # execute
    plot( x , y , properties )

    # show current script
    #BASH("more .temporary.gnu")







#
