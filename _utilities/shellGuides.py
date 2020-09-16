#! /usr/bin/env python





# --- core libraries --------



# - bash
from subprocess import check_output as bashCheck
from subprocess import call as bashCall


# - system reading
import sys





# --- bash interpereting tools -------------------------------------------------



# - acquire arguments
def pullWorkingFile():
    return sys.argv[0]


# - acquire arguments
def pullArguments():
    return sys.argv[1:]


# - execute bash command
def BASH(command):
    if    isinstance(command,list): return str( bashCheck(command) )[2:-3]
    elif  isinstance(command,str):  return str( bashCheck(command.split()) )[2:-3]
    else: print( "could not interperet bash command: {0}".format(command) )

# - execute bash command
def BASHout(command):
    if    isinstance(command,list): return str( bashCall(command) )[2:-3]
    elif  isinstance(command,str):  return str( bashCall(command.split()) )[2:-3]
    else: print( "could not interperet bash command: {0}".format(command) )


# - create escape color string
def bashColor( colorName , *bold ):
# input:
# colorName - "name of text color" string
# bold - "bold text?" 0 or 1 int (optional)
    # manage 'bold' input
    try:
        bold = int(bold[0])
        if bold>0: bold = 1
    except: bold = 0
    # create bash color string dictionary
    colorDictionary = {}
    colorDictionary["BLACK"]   = [ "\033[0;30m" , "\033[1;30m" , "k" , "K" , "black"  , "Black"  , "BLACK"  ]
    colorDictionary["RED"]     = [ "\033[0;31m" , "\033[1;31m" , "r" , "R" , "red"    , "Red"    , "RED"    ]
    colorDictionary["GREEN"]   = [ "\033[0;32m" , "\033[1;32m" , "g" , "G" , "green"  , "Green"  , "GREEN"  ]
    colorDictionary["ORANGE"]  = [ "\033[0;33m" , "\033[1;33m" , "o" , "O" , "orange" , "Orange" , "ORANGE" ]
    colorDictionary["PURPLE"]  = [ "\033[0;34m" , "\033[1;34m" , "m" , "M" , "purple" , "Purple" , "PURPLE" ]
    colorDictionary["PINK"]    = [ "\033[0;35m" , "\033[1;35m" , "p" , "P" , "pink"   , "Pink"   , "PINK"   ]
    colorDictionary["CYAN"]    = [ "\033[0;36m" , "\033[1;36m" , "c" , "C" , "cyan"   , "Cyan"   , "CYAN"   ]
    colorDictionary["GRAY"]    = [ "\033[0;37m" , "\033[1;37m" , "g" , "G" , "gray"   , "grey"   , "Gray" , "Grey" , "GRAY" , "GREY" ]
    # cross-reference 'colorName' with
    for key,colors in colorDictionary.items():
        if colorName in colors:
            return colors[bold]
    # otherwise return [No Color] reset code
    return "\033[0m"





# --- argument interperetation -------------------------------------------------



# - peel off file input
def identifyFileName( variables ):
    # if there are no variables, return current directory
    try:
        if len(variables)==0: return "."
    except:
        return "."
    # iterate through input variables
    for name in variables:
        try:
            # if found, return name
            tmp = open(name,'r')
            tmp.close()
            return name
        except:
            # if permutation of name found, return that
            for j in range(100):
                try:
                    NAME = name + "-" + str(j)
                    tmp = open(NAME,'r')
                    tmp.close()
                    warning("original file '{0}' not found, yet have identified and chosen '{1}'".format(name,NAME))
                    return NAME
                except: pass
    # if none found, return current directory
    return "."


# - interperet variable input
def string2List( String ):
    # 1 - try simple value reading
    try: return [ float(String) ]
    except: pass
    # 2 - check for range
    if not String.find(":")<0:
        idx = String.find(":")
        minValue = int(abs(float( String[:idx]   )))
        maxValue = int(abs(float( String[idx+1:] )))
        return list(range( minValue , maxValue+1 ))
    # 3 - otherwise collect values
    idx = 0
    List = []
    for j in range(len(String)):
        if String[j] not in ["-","0","1","2","3","4","5","6","7","8","9"]:
            List.append( int(float(String[idx:j])) )
            idx = j+1
    List.append( float(String[idx:]) )
    return List

# - convert list into integers
def castListInt( List ):    return [ int(list) for list in List ]


# - convert list into nonnegative integers
def castListUint( List ):   return [ int(abs(list)) for list in List ]


# - prints warnings
def warning(*argv):
    print("\033[1;33m--- WARNING : shellGuides ---\033[0;33m")
    for arg in argv:
        print(arg)
    print("\033[0m")





# ------------------------------------------------------------------------------
