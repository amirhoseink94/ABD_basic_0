CC = g++
CFLAGS = -Wall -g
OPGLFLAGS = -lGL -lGLU -lglut -lGLEW -lglfw -lX11 -lXxf86vm -lXrandr -lpthread -lXi -ldl -lXinerama -lXcursor
HEADERS = Particle.h

res: main.o Particle.o Shape.o Utilities.o
	$(CC) $(CFLAGS) -o res main.o Particle.o Shape.o Utilities.o $(OPGLFLAGS)
#utilities.o: utilities.cpp
#		$(CC) $(CFLAGS) -c utilities.cpp #$(HEADERS)
main.o: main.cpp
	$(CC) $(CFLAGS) -c main.cpp #$(HEADERS)
Utilities.o: src/Utilities.cpp
			$(CC) $(CFLAGS) -c src/Utilities.cpp #$(HEADERS)
Shape.o: src/Shape.cpp
		$(CC) $(CFLAGS) -c src/Shape.cpp #$(HEADERS)
Particle.o: src/Particle.cpp
	$(CC) $(CFLAGS) -c src/Particle.cpp #$(HEADERS)
clean:
	rm -f main.o Particle.o
