# -*- Mode: makefile -*-
# Makefile for GAlib
# Copyright (c) 1996-2005 Matthew Wall, all rights reserved
#
# If you need to customize the build of galib, you should first modify the
# variables in the makevars file.

include galib247/makevars

CXX = g++ -std=c++11
EXEC=main
OBJ=

GA_INC_DIR= galib247
GA_LIB_DIR= galib247/ga

INC_DIRS= -I$(GA_INC_DIR)
LIB_DIRS= -L$(GA_LIB_DIR)

all:
	make clean; make ex ; make run

ex:
	cd galib247; $(MAKE) lib
	$(CXX) $(CXXFLAGS) $(EXEC).cpp $(OBJ) $(INC_DIRS) $(LIB_DIRS) -o $(EXEC) -lga -lm

clean:
	rm -rf *.o *.png *.dat $(EXEC)

.PHONY: clean rebuild

run: $(EXEC)
	./$(EXEC)

