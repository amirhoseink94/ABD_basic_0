CC = g++
CFLAGS = -Wall -g
OPGLFLAGS = -lGL -lGLU -lglut -lGLEW -lglfw -lX11 -lXxf86vm -lXrandr -lpthread -lXi -ldl -lXinerama -lXcursor
ARMAFLAG = -larmadillo
HEADERS = Particle.h

shape_generator: generator.o Shape.o
	$(CC) $(CFLAGS) -o shape_generator generator.o Shape.o $(OPGLFLAGS) $(ARMAFLAG)
Shape.o: Shape.cpp
		$(CC) $(CFLAGS) -c Shape.cpp #$(HEADERS)
generator.o: generator.cpp
	$(CC) $(CFLAGS) -c generator.cpp
clean:
	rm -f main.o Particle.o
