#!/usr/bin/env bash

# clear debug trap
trap "" DEBUG


# ---- USER argument and parameter handling -----------------

# User input (file name optional)
# STUB
dosDir="doscar_parsed"
baseName=".doscar-"

# --- get ion count for file naming ------------------------

IONS=$(head -n 1 DOSCAR | tr -s ' ' | cut -d ' ' -f 2 )
if [ ${IONS} -lt 10 ]; then
	PAD=1
elif [ ${IONS} -lt 100 ]; then
	PAD=2
elif [ ${IONS} -lt 1000 ]; then
	PAD=3
else
	PAD=4
fi

# ---- MAIN -----------------------------------

# split data into separate files
printf "main split\n"
split -a 3 -p "$(sed '6q;d' DOSCAR)" DOSCAR ${baseName}

# create doscar directory
mkdir -p ${dosDir}

# create header
sed '6q;d' DOSCAR >> "${baseName}aaa"
more ${baseName}aaa | tr -s ' ' | sed -e 's/^[ \t]*//' > "${dosDir}/header.txt"
rm ${baseName}aaa

# create main doscar file
tail -n +2 ${baseName}aab | tr -s ' ' | sed -e 's/^[ \t]*//' > "${dosDir}/doscar.dat"
rm ${baseName}aab

# format files (trim white spaces, kill duplicate spaces)
printf "simplifying pdos files\n"
j=0; for dosFile in ${baseName}*; do
	let j++
	#printf "j = ${j}    dosFile = ${dosFile}\n"
	tail -n +2 ${dosFile} | tr -s ' ' | sed -e 's/^[ \t]*//' > ".p${dosFile}"
	mv ".p${dosFile}" "${dosDir}/pdoscar_ion$(printf "%0${PAD}d" ${j}).dat" 
done
rm ${baseName}*
cd ${dosDir}


# --- HANDLE PROJECTED DOS files -------------------------------

# count number of columns in the DOSCAR data (number excludes the energy column)
COLUMNS=$(head -n 1 "pdoscar_ion$(printf "%0${PAD}d" ${j}).dat" | grep -o -e " " | wc -l)

printf "split pdos files by decompositions\n"
# CASES: break into subfiles based on number of columns
j=0;
if [ ${COLUMNS} -eq 3 ]; then
	echo "nope" # ISPIN=1, l-decomposed ### energy s-DOS p-DOS d-DOS
elif [ ${COLUMNS} -eq 9 ]; then
	echo "nope" # ISPIN=1, lm-decomposed ### energy  s  p_y p_z p_x d_{xy} d_{yz} d_{z2-r2} d_{xz} d_{x2-y2}, ...
elif [ ${COLUMNS} -eq 4 ]; then
	echo "nope" # ISPIN=2, no decomposition ### energy     DOS(up) DOS(dwn)  integrated DOS(up) integrated DOS(dwn)
elif [ ${COLUMNS} -eq 6 ]; then
	echo "nope" # ISPIN=2, l-decomposed
elif [ ${COLUMNS} -eq 18 ]; then
	echo "nope" # ISPIN=2, lm-decomposed
elif [ ${COLUMNS} -eq 12 ]; then
	printf "l-decomposed noncollinear DOSCAR detected...\n" # noncolinear, l-decomposed
	printf "The 4 columns in each 'pDOSCAR...dat' file are 'energy s p d' \n\n" >> "header.txt"
	for pdosFile in "p"*; do
		let j++
		PREFIX="$(printf "%0${PAD}d" ${j}).dat"
		cut -d " " -f 1,2,6,10 ${pdosFile} > "pdoscarT_ion${PREFIX}"
		cut -d " " -f 1,3,7,11 ${pdosFile} > "pdoscarX_ion${PREFIX}"
		cut -d " " -f 1,4,8,12 ${pdosFile} > "pdoscarY_ion${PREFIX}"
		cut -d " " -f 1,5,9,13 ${pdosFile} > "pdoscarZ_ion${PREFIX}"
		mv "pdoscarT_ion${PREFIX}" "${pdosFile}"
	done
elif [ ${COLUMNS} -eq 36 ]; then
	printf "lm-decomposed noncollinear DOSCAR detected...\n" # noncolinear, lm-decomposed
	printf "The 10 columns in each 'pDOSCAR...dat' file are 'energy s py pz px dxy dyz dz2 dxz x2-y2' \n\n" >> "header.txt"
	for pdosFile in "p"*; do
		let j++
		PREFIX="$(printf "%0${PAD}d" ${j}).dat"
		cut -d " " -f 1,2,6,10,14,18,22,26,30,34 ${pdosFile} > "pdoscarT_ion${PREFIX}"
		cut -d " " -f 1,3,7,11,15,19,23,27,31,35 ${pdosFile} > "pdoscarX_ion${PREFIX}"
		cut -d " " -f 1,4,8,12,16,20,24,28,32,36 ${pdosFile} > "pdoscarY_ion${PREFIX}"
		cut -d " " -f 1,5,9,13,17,21,25,29,33,37 ${pdosFile} > "pdoscarZ_ion${PREFIX}"
		mv "pdoscarT_ion${PREFIX}" "${pdosFile}"
	done
else
	echo "Nope" # unknown, just leave everything as is ...
fi


printf "finished splitting\nchecking for a poscar...\n"







# --- if POSCAR available, U\use it to create the ion-named files -------------------------------

if test -f "../POSCAR"; then
	# acquire arrays of the ions and their species
	printf "POSCAR found:\n"
	declare -i NUMBERS
	read -ra NUMBERS <<< $(head -n 7 ../POSCAR | tail -n 1 | tr -s ' ' | sed -e 's/^[ \t]*//')
	read -ra SPECIES <<< $(head -n 6 ../POSCAR | tail -n 1 | tr -s ' ' | sed -e 's/^[ \t]*//')
	# compute ions specified in POSCAR
	j=0;
	for (( i=0; i<${#NUMBERS[@]}; i++ )); do
    	j=$(( ${j} + ${NUMBERS[$i]} ))
    done
    printf "$j ions determined:"
    # only modify file names if POSCAR and DOSCAR are compatible
    if [ $j -eq ${IONS} ]; then
    	# loop by species
    	k=0;
    	for (( i=0; i<${#SPECIES[@]}; i++ )); do
    		# loop over number of this species
    		for (( j=0; j<${NUMBERS[$i]}; j++ )); do
    			PREFIX="$( printf "%0${PAD}d" $((${j}+${k}+1)) ).dat"
    			printf "\n ... ${PREFIX}"
    			for pdosFile in *"${PREFIX}"; do
    				printf "\033[0;33m ${pdosFile} \n\033[0m"
    				mv "${pdosFile}" "$( printf ${pdosFile} | cut -d "." -f 1 )_${SPECIES[$i]}.dat"
    			done
    		done
    		k=$(( ${k} + ${j} ))
    	done
    else
    	printf "The DOSCAR and POSCAR in this directory are not compatible!\nThe POSCAR thinks there are ${j} ions, where the DOSCAR thinks there are ${IONS} \n"
    fi
else
	printf "no POSCAR file found :(\n"
fi





# --- if POSCAR available, use it to create the ion totals -------------------------------

LEN=$( head -n 6 ../DOSCAR | tail -n 1 | cut -d "." -f 3- | cut -d " " -f 2 )

#if test -f "../POSCAR"; then
#	printf "this part could take a little while"
#	# loop over species
#	for (( i=0; i<${#SPECIES[@]}; i++ )); do
#		# loop over pdos file lines
#		j=0
#		while [ $j -lt $LEN ]; do
#			let j++
#			# loop over number of this species
 #   		for (( k=0; k<${NUMBERS[$i]}; k++ )); do
  #  			
   # 		done
#		done
#	done
#fi





# restore DEBUG trap
function setWhite { printf "\033[0m" }
trap setWhite DEBUG


# close out
unset IONS
unset PAD
unset COLUMNS
unset dosFile
unset pdosFile
unset PREFIX
unset baseName
unset dosDir
unset SPECIES
unset NUMBERS
cd -
