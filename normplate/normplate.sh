#!/bin/bash

./normplate "$1" | ./cvpersp "$1" "$2"
