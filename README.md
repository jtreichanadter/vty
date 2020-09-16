# vty
(Short for _VASP-Tyler_)

Hello friends!
You are looking at an organized culmination of scripts I've made while parsing and processing VASP data. My hope is to make these open for public use, So ideally they should be as modular and intuitive (user friendly) as I can reasonably muster.

This library is up on github for sharing purposes, and you might want to configure .bashrc to automatically update the directory at the start of each session.



## Current scripts:

### parseDOSCAR.sh (not complete!)
Converts DOSCAR data into a series of pDOS file in the directory 'doscar_parsed' with meta data in 'header.txt', and the main doscar data in 'doscar.dat'. The main file will be broken up based on the type of VASP run it came from (l-decomposed or lm, ispin-1 , etc.). If a POSCAR is found, then pdos filenames will include the ion on the suffix.
**NOT COMPLETE**
Currently this script is just set to work for noncollinear calculations.
It won't fail otherwise, but might not name files in the most intuitive way.

### parseWAVECAR.f90 (currently noncollinear only!)
**Compiling methods**
For gfortran just use _gfortran -o <executable> parseWAVECAR.f90_
For ifort use _ifort -assume byterecl -o <executable> parseWAVECAR.f90_
**Description**
Basically use this to extract a realspace eigenstate from the Kohn-Sham eigenpairs written in the WAVECAR file.
See the file's header for more info on how to use this executable!


