#*	  _________________________________________________________
#*   /                                                        / 
#*  / Smooth Particle hydrodynamics algorithm implementation /
#* /________________________________________________________/
#*/
#*  Copyright Bruno Faria
#*  University of Aveiro 
#*  Department of physics
#*  21/05/2014
#*
# escolhe compilador
CC=icc			
LD=icc
#AR=xiar
# flags a passar ao compilador 
# -O3 -mssse3 -align -xssse3 -axssse3
CFLAGS= -c -Wall -g  $(USER_DEF) -O3  -qopenmp  -xHOST -unroll -funroll-loops -align -parallel  -xavx  #-qopt-report=4 -vec-report=4


# flags de biblioteca
LDFLAGS =  -mkl -lGL  -lGLU -lglut -lGLEW 
# (ficheiros a compilar) para já fica assim!!!
SOURCES = dynamics.c main.c memory.c shaders.c camera.c renderer.c

# extenções dos ficheiros seja objectos .o ou ficheiros de código 
OBJECTS=$(SOURCES:.c=.o)

# nome do ficheiro executavel desejado...
EXECUTABLE=sph

# clean
RM = rm -f

# especifica um target, neste caso todos...
all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC)  $(OBJECTS) -o $@ $(LDFLAGS)

.c.o:
	$(CC)  $< -o $@ $(CFLAGS)

clean: 
	$(RM) *.o

