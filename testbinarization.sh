#!/bin/sh

normplate=bin/normplate;
edfiles=files/hcspecnum/edgedetect;

for file in $edfiles/1/*.png
do
	mkdir $edfiles/res/`basename $file`;
	$normplate $file $edfiles/res/`basename $file`/res.png;

	a="$edfiles/res/`basename $file`";

	cp bot.png "$a/bot.png";
	cp top.png "$a/top.png";
	cp left.png "$a/left.png";
	cp right.png "$a/right.png";
done
