# Makefile for DPS lake problem
CXX = mpiCC
CXXFLAGS = -c -Wall
#INC = -I /apps/modulefiles/standard/core/boost/1.54.0/include
LIBS = -lm
UNAME_S = $(shell uname -s)

ifneq (, $(findstring SunOS, $(UNAME_S)))
    LIBS += -lnsl -lsocket -lresolv
endif

all: UrbanHeatMPI

UrbanHeatMPI: UrbanHeat.o borgms.o mt19937ar.o
	$(CXX) UrbanHeat.o borgms.o mt19937ar.o -o UrbanHeatMPI

UrbanHeat.o: UHMOEA_ms_ma.cpp ./Borg-1.9/borgms.h
	$(CXX) $(CXXFLAGS) $(INC) UHMOEA_ms_ma.cpp -o UrbanHeat.o

#moeaframework.o: ./../borg/moeaframework.c ./../borg/moeaframework.h
#	$(CXX) $(CXXFLAGS) ./../borg/moeaframework.c

borgms.o: ./Borg-1.9 ./Borg-1.9/borgms.h \
	    ./Borg-1.9/borg.h
	$(CXX) $(CXXFLAGS) ./Borg-1.9/borgms.c $(LIBS)

mt19937ar.o: ./Borg-1.9/mt19937ar.c ./Borg-1.9/mt19937ar.h
	$(CXX) $(CXXFLAGS) ./Borg-1.9/mt19937ar.c $(LIBS)

clean:
	rm -rf *.o
	rm UrbanHeatMPI