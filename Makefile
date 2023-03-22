default: sim
all: sim

CXX = g++
CXXFLAGS = -std=c++20 -Wall -O0 -g -I ./

SIMHDRS = models/Species.h models/Field.h models/Distribution.h models/Grid.h models/Simulation.h models/Species1D1V.h\
		  models/Field1D1V.h models/Species2D3V.h models/Field2D3V.h models/Vector2.h models/Vector3.h
SIMOBJS = models/Species.o models/Field.o models/Distribution.o models/Grid.o models/Simulation.o models/Species1D1V.o\
		  models/Field1D1V.o models/Species2D3V.o models/Field2D3V.o models/Vector2.o models/Vector3.o
SIMLIBS = -lboost_program_options -lblas -llapack


HDRS = $(SIMHDRS) $(VISHDRS)
OBJS = $(SIMOBJS) $(VISOBJS)

%.o : %.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -o $@ -c $<

sim: main.o $(SIMOBJS)
	$(CXX) $(CXXFLAGS) -o sim $^ $(SIMLIBS)

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
	$(CXX) $(CXXFLAGS) -o test $^ $(SIMLIBS)

.PHONY: distributionTests
distributionTests: tests/distributionTests.o $(SIMOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(SIMLIBS)

.PHONY: presetDistributionTests
presetDistributionTests: tests/presetDistributionTests.o $(SIMOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(SIMLIBS)

.PHONY: vectorTests
vectorTests: tests/vectorTests.o $(SIMOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(SIMLIBS)

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
