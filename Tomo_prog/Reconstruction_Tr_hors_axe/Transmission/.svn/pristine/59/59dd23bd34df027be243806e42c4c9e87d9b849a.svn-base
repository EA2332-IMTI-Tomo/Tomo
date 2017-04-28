#!/bin/bash


#GLOS (Graphic Library in Open Source), an ANSI Common Lisp OpenGL subset. 
#Copyright (C) 2001 the GLOS development team (http://glos.sourceforge.net) 

source ~/bin/functions.bash


IMG_SUFFIX="pgm"
TEMPO_MS=18



usage() 
{
    echo "usage: $(basename $0) <src dir> <dst dir> "
    echo "simule une acquisition en cours en faisant apparaitre des hologrammes aussi vite que le microscope en produit. Le but est de nourrir un programme de reconstruction en temps-réel."
    echo "src dir: répertoire où sont stockés les hologrames expérimentaux"
    echo "dst dir: répertoire où ils sont copiés au fur et à mesure (-> à écouter par le programme de reconstruction)"
    exit 1
}





if [ $# -ne 2 ]
    
    then
    usage $#
    exit 1

else
    DIR_SRC=$1
    DIR_DST=$2


    echo "src: $DIR_SRC | dst: $DIR_DST | tempo (ms): $TEMPO_MS"


    TEMPO_S=$(echo "$TEMPO_MS * 0.001" | bc)
    NUM_RECORD=1
    NUM_HOLO=0


    IMAGE_FILES=$(ls $DIR_SRC/*.$IMG_SUFFIX) 
    NB_IMAGE_FILES=$(echo $IMAGE_FILES | wc -w)

    FILENAME_ONE=$(ls $DIR_SRC/*.pgm 2> /dev/null | head -n 1)   
    FILENAME_RADIX=$(basename $FILENAME_ONE)
    SESSION_RADIX=$( echo $FILENAME_RADIX | cut -d '-' -f 1 )

    ITER=0
    
    
    BEGIN_T=$(sleepenh 0.001)
    TIME=$BEGIN_T

    echo "simulation start:"
    date

    while [ $ITER -lt $NB_IMAGE_FILES ]
    do
	NUM_HOLO_DISP=$( echo $(expr $NUM_HOLO + 1) | awk '{printf "%03d\n", $0}' )
	#aussi: printf "%03d\n"

	FILE="$DIR_SRC/$SESSION_RADIX-record$NUM_RECORD-$NUM_HOLO_DISP.$IMG_SUFFIX"

	cp $FILE $DIR_DST

	TIME2=$(sleepenh $TIME $TEMPO_S)
	#sleepenh $TEMPO_S >/dev/null
	TIME=$TIME2
	
	echo -ne "num holo: $NUM_RECORD\033[0K\r"; 



	let ITER=ITER+1

	NUM_RECORD=$(expr $ITER / 4 + 1) 
	let NUM_HOLO=NUM_HOLO+1
	NUM_HOLO=$( expr $NUM_HOLO % 4 )


    done

    END_T=$(sleepenh 0.001)
    TIME_TOTAL=$(echo "$END_T - $BEGIN_T" | bc)

    echo "durée simu: $TIME_TOTAL"

    
    echo "simulation stop:"
    date
    
fi
