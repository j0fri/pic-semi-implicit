default: sim
all: sim

CXX = g++
CXXFLAGS = -std=c++20 -Wall -O0 -g -I ./

SIMHDRS = models/Species.h models/Field.h models/Distribution.h models/Grid.h models/Simulation.h models/Species1D1V.h\
		  models/Field1D1V.h
SIMOBJS = models/Species.o models/Field.o models/Distribution.o models/Grid.o models/Simulation.o models/Species1D1V.o\
		  models/Field1D1V.o
SIMLIBS = -lboost_program_options -lblas -llapack


HDRS = $(SIMHDRS) $(VISHDRS)
OBJS = $(SIMOBJS) $(VISOBJS)

%.o : %.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -o $@ -c $<

sim: main.o $(SIMOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(SIMLIBS)

vis: $(VISOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(VISLIBS)

.PHONY: clean
clean:
	-rm -f *.o **/*.o sim vis *Test*

.PHONY: cleanOutput
cleanOutput:
	-rm -f outputs/*

.PHONY: profiler
profiler: sim
	-module load dev-studio
	-rm -r profiling.er
	-collect -o profiling.er ./sim
	-analyzer profiling.er
	
.PHONY: helperTests
helperTests: tests/helperTests.o
	$(CXX) $(CXXFLAGS) -o test $^

.PHONY: distributionTests
distributionTests: tests/distributionTests.o $(SIMOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

.PHONY: presetDistributionTests
presetDistributionTests: tests/presetDistributionTests.o $(SIMOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

.PHONY: landauFile
landauFile: sim
	./sim  --mode 1
	
.PHONY: landau
landau: sim
	./sim  --mode 2
	
.PHONY: twoStreamFile
twoStreamFile: sim
	./sim  --mode 3

.PHONY: twoStream
twoStream: sim
	./sim  --mode 4
