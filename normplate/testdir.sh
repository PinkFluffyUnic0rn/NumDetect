#!/bin/bash

for file in $1/*.png
do
	./normplate "$file" $2/$(basename "$file")
done;
