CC = g++
CFLAGS = -Wall -g
OPGLFLAGS = -lGL -lGLU -lglut -lGLEW -lglfw -lX11 -lXxf86vm -lXrandr -lpthread -lXi -ldl -lXinerama -lXcursor
 
res: main.o Particle.o
	$(CC) $(CFLAGS) -o res main.o Particle.o $(OPGLFLAGS) 
main.o: main.cpp
	echo "this is going to happen 2"
	$(CC) $(CFLAGS) -c main.cpp 
Particle.o: Particle.cpp
	$(CC) $(CFLAGS) -c Particle.cpp 
clean:
	rm -f res main.o Particle.o 
 
