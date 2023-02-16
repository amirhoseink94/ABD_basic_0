CC = g++
CFLAGS = -Wall -g
OPGLFLAGS = -lGL -lGLU -lglut -lGLEW -lglfw -lX11 -lXxf86vm -lXrandr -lpthread -lXi -ldl -lXinerama -lXcursor
ARMAFLAG = -larmadillo
HEADERS = Particle.h

res: main.o Particle.o Segment.o Face.o Body.o Utilities.o
	$(CC) $(CFLAGS) -o res main.o Particle.o Segment.o Face.o Body.o Utilities.o $(OPGLFLAGS) $(ARMAFLAG)
#utilities.o: utilities.cpp
#		$(CC) $(CFLAGS) -c utilities.cpp #$(HEADERS)
main.o: main.cpp
	$(CC) $(CFLAGS) -c main.cpp
Utilities.o: src/Utilities.cpp
	$(CC) $(CFLAGS) -c src/Utilities.cpp
Body.o: src/Body.cpp
	$(CC) $(CFLAGS) -c src/Body.cpp
Face.o: src/Face.cpp
	$(CC) $(CFLAGS) -c src/Face.cpp
Segment.o: src/Segment.cpp
	$(CC) $(CFLAGS) -c src/Segment.cpp
Particle.o: src/Particle.cpp
	$(CC) $(CFLAGS) -c src/Particle.cpp
clean:
	rm -f main.o Particle.o
