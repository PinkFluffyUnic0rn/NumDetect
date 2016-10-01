#!/bin/bash

for file in $2/*.png
do
	filen=1
	echo "processing: $file"

	./numdetect.sh "$1" "$file" res

	for rs in res/*.png
	do
		fname=$(basename $file)
		cp "$rs" "$3/${fname::-4}_$filen.png"
		((++filen))
	done

	rm res/*.png
done
