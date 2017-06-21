CC      = g++
SOURCEDIR = ./
INSTALLDIR = ~/anaconda2/bin
bruteclust: $(SOURCEDIR)/make_clust_brute.cpp
	$(CC) -o bruteclust $(SOURCEDIR)/make_clust_brute.cpp; mv bruteclust $(INSTALLDIR)


kdhsa102: $(SOURCEDIR)/kdhsa102_main.cpp $(SOURCEDIR)/kdhsa102.h
	$(CC) -o kdhsa102 $(SOURCEDIR)/kdhsa102.cpp $(SOURCEDIR)/kdhsa102_main.cpp; mv kdhsa102 $(INSTALLDIR)

6dimprobable: $(SOURCEDIR)/6dim_main.cpp $(SOURCEDIR)/6dimprobable.h
	$(CC) -o 6dimprobable $(SOURCEDIR)/6dimprobable.cpp $(SOURCEDIR)/6dim_main.cpp; mv 6dimprobable $(INSTALLDIR)

all: bruteclust kdhsa102 6dimprobable

clean:
	rm -f bruteclust
	rm -f kdhsa102
	rm -f 6dimprobable
