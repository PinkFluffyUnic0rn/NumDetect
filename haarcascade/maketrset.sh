#!/bin/bash

echo "$(ls -f "$1" | grep .png | wc -l) \
$(ls -f "$2" | grep .png | wc -l)"

$(dirname $0)/printpaths.sh "$1"
$(dirname $0)/printpaths.sh "$2"
