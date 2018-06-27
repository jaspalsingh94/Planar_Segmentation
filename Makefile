########################################
##
## Makefile
## LINUX compilation
##
##############################################





#FLAGS
C++FLAG = -g -std=c++11 -Wall
MATH_LIBS = -lm

EXEC_DIR=.

.cc.o:
	g++ $(C++FLAG) $(INCLUDES)  -c $< -o $@

#Including
INCLUDES=  -I.

#-->All libraries (without LEDA)
LIBS_ALL =  -L/usr/lib -L/usr/local/lib


Cpp_OBJ1= main.cpp DisjSets.o
PROGRAM_1=p1
$(PROGRAM_1): $(Cpp_OBJ1)
	g++ $(C++FLAG) -o $(EXEC_DIR)/$@ $(Cpp_OBJ1) $(INCLUDES) $(LIBS_ALL)

all:
	make $(PROGRAM_1)


clean:
	(rm -f *.o;)

(:
