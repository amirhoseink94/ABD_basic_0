CC = g++
CFLAGS = -Wall -g
OPGLFLAGS = -lGL -lGLU -lglut -lGLEW -lglfw -lX11 -lXxf86vm -lXrandr -lpthread -lXi -ldl -lXinerama -lXcursor
ARMAFLAG = -larmadillo
extFLAG = -I ~/Desktop/ThreeDEnvironment/ext/
HEADERS = Particle.h

res: main.o CCD
	$(CC) $(CFLAGS) $(extFLAG) -o res *.o #main.o Particle.o Segment.o Face.o Body.o Engin.o ccd.cpp.o avx.cpp.o interval_root_finder.cpp.o interval.cpp.o
#CCD:
#		$(CC) $(CFLAGS) -c ext/tight_inclusion/*.cpp #$(HEADERS)
main.o: main.cpp
	$(CC) $(extFLAG) $(CFLAGS) -c main.cpp
CCD: 
	$(CC) $(CFLAGS) $(extFLAG) -c ccd.cpp
	$(CC) $(CFLAGS) $(extFLAG) -c avx.cpp
	$(CC) $(CFLAGS) $(extFLAG) -c interval.cpp
	$(CC) $(CFLAGS) $(extFLAG) -c interval_root_finder.cpp

