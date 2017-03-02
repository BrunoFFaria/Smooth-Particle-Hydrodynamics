#*	  ________________________________________________
#*   /                                               / 
#*  / Smooth particle hydrodynamics implementation  /
#* /_______________________________________________/
#*/
#*  Copyright Bruno Faria
#*  University of Aveiro 
#*  Department of physics
#*  21/01/2017
#*
# choose compiler (I use ICC)
CC=icc			
LD=icc
#AR=xiar
# flags to pass the compiler
# -O3 -mssse3 -align -xssse3 -axssse3
CFLAGS= -c -Wall -g  $(USER_DEF) -O3 -xHOST -qopenmp -unroll -funroll-loops -align -parallel  -xavx  #-qopt-report=4 -vec-report=4


# flags of libraries
LDFLAGS =  -mkl -lGL  -lGLU -lglut -lGLEW 
# (files to compile)
SOURCES = dynamics.c main.c marching_cubes.c memory.c shaders.c zpr.c

# rule for getting an object .o 
OBJECTS=$(SOURCES:.c=.o)

# target name
EXECUTABLE=sph

# clean
RM = rm -f

# 
all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC)  $(OBJECTS) -o $@ $(LDFLAGS)

.c.o:
	$(CC)  $< -o $@ $(CFLAGS)

clean: 
	$(RM) *.o

