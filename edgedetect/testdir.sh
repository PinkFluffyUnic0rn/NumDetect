#!/bin/bash

for file in $1/*.png
do
	$(dirname $0)/../bin/normplate "$file" $2/$(basename "$file")
done;
