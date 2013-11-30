CC = g++
ifeq ($(shell sw_vers 2>/dev/null | grep Mac | awk '{ print $$2}'),Mac)
	CFLAGS = -g -fopenmp -DGL_GLEXT_PROTOTYPES -I./include/ -I/usr/X11/include -DOSX -O2
	LDFLAGS = -framework GLUT -framework OpenGL \
    	-L"/System/Library/Frameworks/OpenGL.framework/Libraries" \
    	-lGL -lGLU -lm -lstdc++
else
	CFLAGS = -msse4 -g -fopenmp -DGL_GLEXT_PROTOTYPES -Iglut-3.7.6-bin 
	LDFLAGS = -msse4 -fopenmp -lglut -lGLU
endif
	
RM = /bin/rm -f 
all: main 
main: primitives.o shapes.o loadscene.o scene.o lodepng.o
	$(CC) $(CFLAGS) -o scene primitives.o shapes.o scene.o loadscene.o lodepng.o  $(LDFLAGS) 
scene.o: scene.cpp
	$(CC) $(CFLAGS) -c scene.cpp -o scene.o
loadscene.o: src/loadscene.cpp
	$(CC) $(CFLAGS) -c src/loadscene.cpp -o loadscene.o
primitives.o: src/primitives.cpp
	$(CC) $(CFLAGS) -c src/primitives.cpp -o primitives.o
shapes.o: src/shapes.cpp
	$(CC) $(CFLAGS) -c src/shapes.cpp -o shapes.o
lodepng.o: src/imgwriter/lodepng.cpp
	$(CC) $(CFLAGS) -c src/imgwriter/lodepng.cpp -o lodepng.o
clean: 
	$(RM) *.o test scene
 
 


