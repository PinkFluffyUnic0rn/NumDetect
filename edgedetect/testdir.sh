#!/bin/bash

for file in $1/*.png
do
	$(dirname $0)/../bin/normplate "$file" "$2" "$3" $4/$(basename "$file")
done;
