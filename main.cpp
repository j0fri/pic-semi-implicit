#include <boost/program_options.hpp>

#include "models/Config.h"
#include "models/Simulation.h"
#include "helpers/preset_configs.h"

//namespace po = boost::program_options;

int main(int argc, char* argv[]){
//	po::options_description opts("Available options.");
//	opts.add_options()
//		("mode", po::value<int>(), "Input preset mode.")
//		("saveParticles", po::value<bool>()->default_value(false), "Save particle position.")
//		("saveEnergy", po::value<bool>()->default_value(true), "Save field and kinetic energy.")
//		("saveInitial", po::value<bool>()->default_value(false), "Save initial electric field and charge distribution.")
//		("saveInitialVelocities", po::value<bool>()->default_value(false), "Save initial velocity distribution.")
//		("loadInitialParticles", po::value<bool>()->default_value(false), "Load initial positions and velocities from file.")
//		("loadInitialField", po::value<bool>()->default_value(false), "Load initial field from file.")
//		("wavePertMag", po::value<float>()->default_value(0.0), "Magnitude of wave perturbation in electric field.")
//		("ePertMag", po::value<float>()->default_value(0.0), "Magnitude of perturbation in distribution of electrons.")
//		("Nx", po::value<int>()->default_value(30), "Field grid discretization points in x.")
//		("Np", po::value<int>()->default_value(10000), "Number of particles.")
//		("Lx", po::value<float>()->default_value(1), "Field non-dimensional length in x.")
//		("Kb", po::value<float>()->default_value(0.00001), "Boltzmann constant.")
//		("T0", po::value<float>()->default_value(1.0), "Initial temperature (K).")
//		("C", po::value<float>()->default_value(1.0), "Speed of light.")
//		("e0", po::value<float>()->default_value(1.0), "Permittivity.")
//		("T", po::value<float>()->default_value(10.0), "Total time.")
//		("dt", po::value<float>()->default_value(0.01), "Time step.")
//		("qe", po::value<float>()->default_value(-1.0), "Electron charge.")
//		("qi", po::value<float>()->default_value(1.0), "Ion charge.")
//		("me", po::value<float>()->default_value(1.0), "Electron mass")
//		("mi", po::value<float>()->default_value(2000.0), "Ion mass")
//		("ue", po::value<float>()->default_value(0.0), "Electron average velocity")
//		("ui", po::value<float>()->default_value(0.0), "Ion average velocity")
//		("rho0", po::value<float>()->default_value(1.0), "Initial unperturbed charge density");
//
//    po::variables_map vm;
//    po::store(po::parse_command_line(argc, argv, opts), vm);

    //auto config = preset_configs::constVelocityX<double>();
//    auto config = preset_configs::landau2D3VX<double>(60,3);
    //auto config = preset_configs::landau2D3VXWave<double>(30,2);
    //auto config = preset_configs::landau2D3VXWaveStatic<double>(30,2);
    //auto config = preset_configs::magneticGyration<double>();
    //auto config = preset_configs::magneticGyrationX<double>();
    //auto config = preset_configs::constPotentialWell<double>();
    //auto config = preset_configs::constPotentialWellFile<double>();
    //auto config = preset_configs::electronBeam<double>(50,50);
    //auto config = preset_configs::fiveParticles<double>();
    //auto config = preset_configs::diode<double>(100000,1,2,10,10);
    auto config = preset_configs::langmuir(1000000,60,3,0.1);
//    config.saveConfig.saveSolverSteps = true;
//    config.fieldConfig.solverTolerance = 1e-2;
    config.saveConfig=preset_save_configs::None<double,2,3>();

	Simulation<double,2,3> sim(config);
    sim.initialise();
    sim.run();
}