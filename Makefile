default: sim
all: sim

CXX = g++
CXXFLAGS = -std=c++17 -Wall -O3 -g -I ./

SIMHDRS = models/Species.h models/Field.h models/Distribution.h models/Grid.h models/Simulation.h models/Species1D1V.h\
		  models/Field1D1V.h models/Species2D3V.h models/Field2D3V.h models/Vector2.h models/Vector3.h\
		  models/DistributionGrid.h

SIMOBJS = models/Species.o models/Field.o models/Distribution.o models/Grid.o models/Simulation.o models/Species1D1V.o\
		  models/Field1D1V.o models/Species2D3V.o models/Field2D3V.o models/Vector2.o models/Vector3.o\
		  models/DistributionGrid.cpp

SIMLIBS = -lboost_program_options -lblas -llapack

HELPERHDRS = helpers/math_helper.h helpers/string_helper.h helpers/output_helper.h

HELPERCPPS = helpers/math_helper.cpp helpers/string_helper.cpp helpers/output_helper.cpp

PRESETHDRS = helpers/preset_configs.h helpers/preset_distributions.h helpers/preset_fields.h helpers/preset_species.h

PRESETCPPS = helpers/preset_configs.cpp helpers/preset_distributions.cpp helpers/preset_fields.cpp helpers/preset_species.cpp

TESTHDRS = testModels/Field2D3VConst.h

TESTOBJS = testModels/Field2D3VConst.o

HDRS = $(SIMHDRS) $(VISHDRS) $(TESTHDRS) $(HELPERHDRS)

%.o : %.cpp $(HDRS) $(HELPERCPPS)
	$(CXX) $(CXXFLAGS) -o $@ -c $<

main.o: main.cpp $(HDRS) $(PRESETHDRS) $(HELPERCPPS) $(PRESETHDRS)

sim: main.o $(SIMOBJS)
	$(CXX) $(CXXFLAGS) -o sim $^ $(SIMLIBS)

.PHONY: clean
clean:
	-rm -f *.o **/*.o sim vis *Test* test

.PHONY: cleanOutput
cleanOutput:
	-rm -f outputs/*

#.PHONY: profiler
#profiler: sim
#	-module load dev-studio
#	-rm -r profiling.er
#	-collect -o profiling.er ./sim
#	-analyzer profiling.er
	
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
vectorTests: tests/vectorTests.o models/Vector2.o models/Vector3.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(SIMLIBS)

.PHONY: constFieldTest1
constFieldTest1: tests/constFieldTest1.o $(SIMOBJS) $(TESTOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(SIMLIBS)

.PHONY: samePotentialWellTest
samePotentialWellTest: tests/samePotentialWellTest.o $(SIMOBJS) $(TESTOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(SIMLIBS)

.PHONY: sameLandauTest
sameLandauTest: tests/sameLandauTest.o $(SIMOBJS) $(TESTOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(SIMLIBS)

export PATH := /home/jf1519/Downloads/tmp/OracleDeveloperStudio12.5-linux-x86-bin/developerstudio12.5/bin:$(PATH)
export MANPATH := /home/jf1519/Downloads/tmp/OracleDeveloperStudio12.5-linux-x86-bin/developerstudio12.5/man:$(MANPATH)

.PHONY: profiler
profiler: sim
	-rm -r profiler.er
	collect -o profiler.er ./sim
	analyzer profiler.er


