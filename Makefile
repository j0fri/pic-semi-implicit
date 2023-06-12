default: sim
all: sim

CXX = mpic++
CXXFLAGS = -std=c++17 -Wall -O3 -g -I ./

SIMHDRS = models/Species.h models/Field.h models/Distribution.h models/Grid.h models/Simulation.h models/Species1D1V.h\
		  models/Field1D1V.h models/Species2D3V.h models/Field2D3V.h models/Vector2.h models/Vector3.h\
		  models/DistributionGrid.h models/Field2D3VExplicit.h

SIMOBJS = models/Species.o models/Field.o models/Distribution.o models/Grid.o models/Simulation.o models/Species1D1V.o\
		  models/Field1D1V.o models/Species2D3V.o models/Field2D3V.o models/Vector2.o models/Vector3.o\
		  models/DistributionGrid.o models/Field2D3VExplicit.o

SIMLIBS = -lblas -llapack

HELPERHDRS = helpers/math_helper.h helpers/string_helper.h helpers/output_helper.h

HELPERCPPS = helpers/math_helper.cpp helpers/string_helper.cpp helpers/output_helper.cpp

PRESETHDRS = helpers/preset_configs.h helpers/preset_distributions.h helpers/preset_fields.h helpers/preset_species.h\
			 helpers/preset_save_configs.h

PRESETCPPS = helpers/preset_configs.cpp helpers/preset_distributions.cpp helpers/preset_fields.cpp helpers/preset_species.cpp\
			 helpers/preset_save_configs.h

TESTHDRS = testModels/Field2D3VConst.h

TESTOBJS = testModels/Field2D3VConst.o

HDRS = $(SIMHDRS) $(VISHDRS) $(TESTHDRS) $(HELPERHDRS)

%.o : %.cpp $(HDRS) $(HELPERCPPS)
	$(CXX) $(CXXFLAGS) -o $@ -c $<

main.o: main.cpp $(HDRS) $(PRESETHDRS) $(HELPERCPPS) $(PRESETHDRS) $(PRESETCPPS)

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

.PHONY: benchmark1
benchmark1: tests/benchmark1.o $(SIMOBJS) $(TESTOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(SIMLIBS)

.PHONY: benchmark2
benchmark2: tests/benchmark2.o $(SIMOBJS) $(TESTOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(SIMLIBS)

.PHONY: benchmark3
benchmark3: tests/benchmark3.o $(SIMOBJS) $(TESTOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(SIMLIBS)

.PHONY: langmuir
langmuir: tests/langmuir.o $(SIMOBJS) $(TESTOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(SIMLIBS)

.PHONY: langmuir2
langmuir2: tests/langmuir2.o $(SIMOBJS) $(TESTOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(SIMLIBS)

.PHONY: langmuir3
langmuir3: tests/langmuir3.o $(SIMOBJS) $(TESTOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(SIMLIBS)

.PHONY: langmuir4
langmuir4: tests/langmuir4.o $(SIMOBJS) $(TESTOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(SIMLIBS)

.PHONY: runtimeNp
runtimeNp: tests/runtimeNp.o $(SIMOBJS) $(TESTOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(SIMLIBS)

.PHONY: runtimeNg
runtimeNg: tests/runtimeNg.o $(SIMOBJS) $(TESTOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(SIMLIBS)

.PHONY: hotspots
hotspots: tests/hotspots.o $(SIMOBJS) $(TESTOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(SIMLIBS)

.PHONY: optratio
optratio: tests/optratio.o $(SIMOBJS) $(TESTOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(SIMLIBS)

.PHONY: optratio2
optratio2: tests/optratio2.o $(SIMOBJS) $(TESTOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(SIMLIBS)

.PHONY: optaccuracy
optaccuracy: tests/optaccuracy.o $(SIMOBJS) $(TESTOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(SIMLIBS)

.PHONY: standalone
standalone: main.o $(SIMOBJS)
	$(CXX) $(CXXFLAGS) -static -o sim $^ $(SIMLIBS)

export PATH := /home/jf1519/Downloads/tmp/OracleDeveloperStudio12.5-linux-x86-bin/developerstudio12.5/bin:$(PATH)
export MANPATH := /home/jf1519/Downloads/tmp/OracleDeveloperStudio12.5-linux-x86-bin/developerstudio12.5/man:$(MANPATH)

.PHONY: profiler
profiler: sim
	-rm -r profiler.er
	nice -20 collect -o profiler.er -M OMPT mpiexec --np 2 --bind-to none -- ./sim
	analyzer profiler.er

.PHONY: HPCprofiler
HPCprofiler: sim
	-rm -r profiler.er
	collect -o profiler.er -M OMPT mpiexec --np 8 --bind-to none -- ./sim
	analyzer profiler.er