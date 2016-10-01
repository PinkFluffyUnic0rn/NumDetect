#!/bin/bash

for file in $1/*.png
do
	./normplate.sh "$file" $2/$(basename "$file")
done;
