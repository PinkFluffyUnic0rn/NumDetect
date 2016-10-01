#!/bin/bash

for file in $1/*.png
do
	./segplate "$file" $2/$(basename "$file")
done;
