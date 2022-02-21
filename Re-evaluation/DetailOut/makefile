# Makefile for DPS lake problem
CC = g++
FLAGS = -Wall -O3 -Wno-unused-local-typedefs
INCL = -I /apps/modulefiles/standard/core/boost/1.54.0/include

SOURCES = $(wildcard *.cpp)
OBJECTS = $(SOURCES:.cpp=.o)
EXE = UrbanHeat_RE_detail.exe

all: $(SOURCES) $(EXE)
	rm $(OBJECTS)

.cpp.o:
	$(CC) -c $(CFLAGS) $^ -o $@ $(INCL)
	
$(EXE): $(OBJECTS)
	$(CC) $(OBJECTS) $(CFLAGS) -o $@ $(INCL)

clean:
	rm -f $(OBJECTS) $(EXE)
