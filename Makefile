default: sim

CXX = g++
CXXFLAGS = -Wall -O0 -g

SIMHDRS = Species.h
SIMOBJS = main.o Species.o
SIMLIBS = -lboost_program_options -lblas -llapack

VISHDRS = 
VISOBJS = vis.o
VISLIBS = 

HDRS = $(SIMHDRS) $(VISHDRS)
OBJS = $(SIMOBJS) $(VISOBJS)

%.o : %.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -o $@ -c $<

sim: $(SIMOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(SIMLIBS)

vis: $(VISOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(VISLIBS)

.PHONY: clean
clean:
	-rm -f *.o sim vis

.PHONY: profiler
profiler: sim
	-module load dev-studio
	-rm -r profiling.er
	-collect -o profiling.er ./sim
	-analyzer profiling.er

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
