default: sim

CXX = g++

CXXFLAGS = -std=c++20 -Wall -O0 -g -I ./

SIMHDRS = Species.h Distribution.h helpers/math_helper.h


SIMOBJS = main.o Species.o Distribution.o
SIMLIBS = -lboost_program_options -lblas -llapack
SIMDEPS = math_helper.cpp

VISHDRS = 
VISOBJS = vis.o
VISLIBS = 

HDRS = $(SIMHDRS) $(VISHDRS)
OBJS = $(SIMOBJS) $(VISOBJS)

%.o : %.cpp $(HDRS) $(SIMDEPS)
	$(CXX) $(CXXFLAGS) -o $@ -c $<

sim: $(SIMOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(SIMLIBS)

vis: $(VISOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(VISLIBS)

.PHONY: clean
clean:
	-rm -f **/*.o sim vis test

.PHONY: profiler
profiler: sim
	-module load dev-studio
	-rm -r profiling.er
	-collect -o profiling.er ./sim
	-analyzer profiling.er

.PHONY: distributionTests
distributionTests: tests/distributionTests.o Distribution.o
	$(CXX) $(CXXFLAGS) -o test $^
	
.PHONY: helperTests $(SIMDEPS)
helperTests: tests/helperTests.o
	$(CXX) $(CXXFLAGS) -o test $^

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
