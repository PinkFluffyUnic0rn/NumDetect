#!/bin/sh

normplate=bin/normplate;
edfiles=files/hccarnum/edgedetect;

for file in $edfiles/nums/*.png
do
	mkdir $edfiles/res/`basename $file`;
	$normplate $file edgedetect/model.png \
	$edfiles/res/`basename $file`/res.png;

	a="$edfiles/res/`basename $file`";

	cp bot.png "$a/bot.png";
	cp top.png "$a/top.png";
	cp left.png "$a/left.png";
	cp right.png "$a/right.png";
done
