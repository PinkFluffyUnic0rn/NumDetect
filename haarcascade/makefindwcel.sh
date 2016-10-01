#!/bin/bash

echo "$(ls -f "$1" | grep .png | wc -l) \
$(ls -f "$2" | grep .png | wc -l)"

./printpaths.sh "$1"
./printpaths.sh "$2"
