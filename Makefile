DIRS = edgedetect haarcascade neuronet readvideo segplate
DIRSMAKE = $(DIRS:%=make-%)
DIRSCLEAN = $(DIRS:%=clean-%)

all : $(DIRSMAKE)

.PHONY : $(DIRSMAKE)

$(DIRSMAKE) :
	cd $(@:make-%=%); \
	make; \
	cd ..\

.PHONY : $(DIRSCLEAN)

$(DIRSCLEAN) :
	cd $(@:clean-%=%); \
	make clean; \
	cd ..\

.PHONY : clean

clean : $(DIRSCLEAN)
