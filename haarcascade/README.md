# HaarCascade
Haar cascade:
------

Initilizing Haar cascade:
	./initcascade [window width] [window height] [features directory] [output file]

Finding weak classifiers:
	./maketrset.sh [good examples directory] [bad examples directory]
	./findwcs [trainig set] [initilized cascade, that can already contain wcs] [output file] [number or wcs to find]

Building cascade:

	./maketrset.sh [test examples directory] [bad examples directory]
	./buildcascade [training set] [cascade with found wcs] [output file]

Scaning image with Haar cascade:

	./scanimage [cascade] [input image] [output directory]
