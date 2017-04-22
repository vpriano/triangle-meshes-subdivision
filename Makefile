## Computer Graphics
## Makefile
##
#########################################

LIBS = -lglut -lGLU -lGL -lXmu -lXext -lXi -lX11 -lm
CC = g++

## Object files and executables
MAIN_OUT = a.out

## Requirements for each command
MAIN_REQS = main.cpp

## Targets to compile for each command
MAIN_TARGETS = main.cpp

all: main

## Main 
main: $(MAIN_REQS)
        $(CC) $(MAIN_TARGETS) $(LIBS) -o $(MAIN_OUT)

clean:
        rm -f *~ *.o *.out
