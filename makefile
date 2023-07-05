CC = g++
CFLAGS_NW = -g
CFLAGS = -Wall -g
OPGLFLAGS = -lGL -lGLU -lglut -lGLEW -lglfw -lX11 -lXxf86vm -lXrandr -lpthread -lXi -ldl -lXinerama -lXcursor
ARMAFLAG = -larmadillo
AVCFLAG = -lavcodec -lavformat -lavutil -lswscale
extFLAG = -I ~/Desktop/ThreeDEnvironment/ext/
HEADERS = Particle.h
CVFLAGS = `pkg-config --cflags --libs opencv4`

res: main.o Particle.o Segment.o Face.o Body.o Engin.o ccd.o avx.o internal.o internal_root_finder.o
	$(CC) $(CFLAGS) $(extFLAG) -o res *.o $(OPGLFLAGS) $(ARMAFLAG) $(CVFLAGS)
main.o: main.cpp
	$(CC) $(extFLAG) $(CFLAGS) $(CVFLAGS) -c main.cpp
Engin.o: src/Engin.cpp
	$(CC) $(CFLAGS) $(extFLAG) -c src/Engin.cpp
#Utilities.o: src/Utilities.cpp
#	$(CC) $(CFLAGS) -c src/Utilities.cpp
Body.o: src/Body.cpp
	$(CC) $(CFLAGS) -c src/Body.cpp
Face.o: src/Face.cpp
	$(CC) $(CFLAGS) -c src/Face.cpp
Segment.o: src/Segment.cpp
	$(CC) $(CFLAGS) -c src/Segment.cpp
Particle.o: src/Particle.cpp
	$(CC) $(CFLAGS) -c src/Particle.cpp
ccd.o: CCD/ccd.cpp
	$(CC) $(CFLAGS_NW) $(extFLAG) -c CCD/ccd.cpp
avx.o: CCD/avx.cpp
	$(CC) $(CFLAGS_NW) $(extFLAG) -c CCD/avx.cpp
internal.o: CCD/interval.cpp
	$(CC) $(CFLAGS_NW) $(extFLAG) -c CCD/interval.cpp
internal_root_finder.o: CCD/interval_root_finder.cpp
	$(CC) $(CFLAGS_NW) $(extFLAG) -c CCD/interval_root_finder.cpp
clean:
	rm -f main.o Particle.o
