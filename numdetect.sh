#!/bin/bash

./scanimage $1 $2 $3

if [ $? -ne 0 ]
then
	exit $?
fi

#for plate in $3/*.png
#do
#	./normplate.sh "$plate" "$plate"
	
#	convert -resize x20 "$plate" "$plate"

#	./numseg seg.hc "$plate" tmp

#	if [ $? -ne 0 ]
#	then
#		exit $?
#	fi

#	for sym in tmp/*.png
#	do
#		./readsym "nn" "$sym"
#	done
	
#	printf "\n"
#done

#rm $3/*.png
#rm tmp/*.png
