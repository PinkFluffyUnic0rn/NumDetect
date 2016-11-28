#!/bin/bash

printf "%d %d\n" "9" "20"
printf "%d %d %d\n" "2" "25" "$(find $1/* -type d | wc -l)"

printf "%d %d %d\n" "$(find $1/*/*.png -type f | wc -l)" "$(find $1/* -type d | wc -l)" "$2"

dname=$(find $1/* -type d -exec basename {} \; | sort -k1)

dirn=0
for dir in $dname
do
	printf "%s\n" "-$dir"

	for file in $1/$dir/*.png
	do
		printf "%s\n" "$file"
	done

	((++dirn))
done
