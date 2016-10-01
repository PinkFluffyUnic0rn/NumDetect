#!/bin/bash

filen=0
for file in $1/*.png
do 
	((++filen))
	
	echo "$file"
done
