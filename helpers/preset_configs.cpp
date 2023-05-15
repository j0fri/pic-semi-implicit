template <typename T>
Config<T,1,1> preset_configs::landau1D1V(){
     Config<T,1,1> config{
        std::vector<typename Config<T,1,1>::SpeciesConfig>{{
               preset_species::Uniform1D1V<T>(50000,1,-1,4*M_PI,30,1,1),
               preset_species::Uniform1D1V<T>(50000,2000,1,4*M_PI,30,1,1),
        }},
        preset_fields::Default1D1V<T>(4*M_PI,30,1,1),
        {10,0.01},
        {false, true, false, false, true, false, false, true, false, 0.1, "outputs/"},
        {{true},{}}, //Only periodic boundary conditions
        true,
        false
    };
    //Add perturbation in electrons
    config.speciesConfig[0].xDist += preset_distributions::Sin<T,1>(0,0.2,(T)1.0/2,0);
    return config;
}

template <typename T>
Config<T,2,3> preset_configs::constVelocityX(){
    Config<T,2,3> config{
            std::vector<typename Config<T,2,3>::SpeciesConfig>{{
            //preset_species::Uniform2D3V<T>(50000,1,-1,1,1,10,10,1,0),
            preset_species::TopHat2D3V<T>(50000,1,-1,0,1,0,1),
            preset_species::TopHat2D3V<T>(50000,2000,1,0,1,0,1),
    }},
    preset_fields::Default2D3V<T>(1,1,10,10,1,1),
    {10,0.1},
    {false, true, false, false, true, true, true, true, true, 0.1, "outputs/"},
    {{true,true},{}}, //Only periodic boundary conditions
    true,
    true
    };
    //Add constant velocity to electrons:
    config.speciesConfig[0].initialVGrid.dimensions[0] = {1,1,1}; //Add constant x-velocity
    config.fieldConfig.initialiseFromSpecies = false;
    return config;
}

template <typename T>
Config<T,2,3> preset_configs::landau2D3VX(unsigned int Nx, unsigned int Ny){
    Config<T,2,3> config{
        std::vector<typename Config<T,2,3>::SpeciesConfig>{{
               preset_species::Uniform2D3V<T>(50000,1,-1,1,1,Nx,Ny,1,0.001),
               preset_species::Uniform2D3V<T>(50000,2000,1,1,1,Nx,Ny,1,0.001),
        }},
        preset_fields::Default2D3V<T>(1,1,Nx,Ny,1,1),
        {10,0.1},
        {false, true, false, false, true, true, true, true, true, 0.1, "outputs/"},
        {{true,true},{}}, //Only periodic boundary conditions
        true,
        true
    };
    //Add perturbation in electrons
    config.speciesConfig[0].xDist += preset_distributions::Sin<T,2>(0,0.2,(T)2*M_PI,0);

    return config;
}

template <typename T>
Config<T,2,3> preset_configs::landauFile(){
    Config<T,2,3> config{
            std::vector<typename Config<T,2,3>::SpeciesConfig>{{
                   preset_species::fromFile<T>(50000,1,-1,"./inputs/landau/electronPositions.txt",
                                                    "./inputs/landau/electronVelocities.txt"),
                   preset_species::fromFile<T>(50000,2000,1,"./inputs/landau/ionPositions.txt",
                                                   "./inputs/landau/ionVelocities.txt")
           }},
            preset_fields::Default2D3V<T>(1,1,30,5,1,1),
            {10,0.1},
            {false, true, false, false, true, true, true, true, true, 0.1, "outputs/"},
            {{true,true},{}}, //Only periodic boundary conditions
            true,
            true
    };
    return config;
}

template <typename T>
Config<T,2,3> preset_configs::landau2D3VXWave(unsigned int Nx, unsigned int Ny){
    Config<T,2,3> config{
            std::vector<typename Config<T,2,3>::SpeciesConfig>{{
                   preset_species::Uniform2D3V<T>(50000,1,-1,1,1,Nx,Ny,1,0.001),
                   preset_species::Uniform2D3V<T>(50000,2000,1,1,1,Nx,Ny,1,0.001),
           }},
            preset_fields::Default2D3V<T>(1,1,Nx,Ny,1,1),
            {10,0.1},
            {false, true, false, false, true, true, true, true, true, 0.1, "outputs/"},
            {{true,true},{}}, //Only periodic boundary conditions
            true,
            true
    };
    //Add perturbation in electric field
    config.fieldConfig.forcedE[0] += preset_distributions::Sin<T,2>(0,1,(T)2*M_PI,0);

    return config;
}

template <typename T>
Config<T,2,3> preset_configs::constAccelerationX(){
    Config<T,2,3> config{
            std::vector<typename Config<T,2,3>::SpeciesConfig>{{
                    preset_species::TopHat2D3V<double>(50000, 1, 1,-0.2, 0.2,-0.2, 0.2)
            }},
            preset_fields::ConstE2D3V<T>(0,(T)1,20,20),
            {10,0.1},
            {false, true, false, false, false, true, true, false, true, 0.1, "outputs/"},
            {{true,true},{}}, //Only periodic boundary conditions
            true,
            true
    };
    return config;
}

template <typename T>
Config<T,2,3> preset_configs::magneticGyration(){
    Config<T,2,3> config{
            std::vector<typename Config<T,2,3>::SpeciesConfig>{{
                   preset_species::TopHat2D3V<double>(50000, 1, 1,-0.2, 0.2,0.3, 0.7)
           }},
            preset_fields::ConstB2D3V<T>(2,(T)1,20,20),
            {10,0.1},
            {false, true, false, false, false, true, true, false, true, 0.1, "outputs/"},
            {{true,true},{}}, //Only periodic boundary conditions
            true,
            true
    };
    config.speciesConfig[0].initialVGrid.dimensions[0] = {0.5,0.5,1}; //Add constant x-velocity
    return config;
}

template <typename T>
Config<T,2,3> preset_configs::magneticGyrationX(){
    Config<T,2,3> config{
            std::vector<typename Config<T,2,3>::SpeciesConfig>{{
                   preset_species::TopHat2D3V<double>(50000, 1, 1,-0.2, 0.2,0.3, 0.7)
            }},
            preset_fields::ConstB2D3V<T>(0,(T)1,20,20),
            {10,0.1},
            {false, true, false, false, false, true, true, false, true, 0.1, "outputs/"},
            {{true,true},{}}, //Only periodic boundary conditions
            true,
            true
    };
    config.speciesConfig[0].initialVGrid.dimensions[2] = {-0.5,-0.5,1}; //Add constant z-velocity
    return config;
}

template <typename T>
Config<T,2,3> preset_configs::constPotentialWell(){
    Config<T,2,3> config{
            std::vector<typename Config<T,2,3>::SpeciesConfig>{{
                   preset_species::TopHat2D3VBoltzmann<double>(50000, 1, 1,-0.2, 0.2,0.3, 0.7,1,0.01)
            }},
            preset_fields::ConstE2D3V<T>(0,(T)0,20,20),
            {100,0.1},
            {false, true, false, false, true, true, true, true, false, 0.1, "outputs/"},
            {{true,true},{}}, //Only periodic boundary conditions
            true,
            true
    };
    config.fieldConfig.forcedE[0].f = ([](const std::array<T,2>& arr){return -arr[0];}); //Ex = -x
    config.fieldConfig.forcedE[1].f = ([](const std::array<T,2>& arr){return -arr[1];}); //Ey = -y
    config.speciesConfig[0].initialVGrid.dimensions[0] = {0.2,0.2,1}; //Add constant x-velocity
    return config;
}

template <typename T>
Config<T,2,3> preset_configs::constPotentialWellFile(){
    Config<T,2,3> config{
            std::vector<typename Config<T,2,3>::SpeciesConfig>{{
                preset_species::fromFile<double>(50000,1,1,"./inputs/constPotentialWell/speciesPosition.txt",
                                                 "./inputs/constPotentialWell/speciesVelocity.txt")
            }},
            preset_fields::ConstE2D3V<T>(0,(T)0,20,20),
            {10,0.1},
            {true, true, true, false, false, true, true, false, true, 10, "outputs/"},
            {{true,true},{}}, //Only periodic boundary conditions
            true,
            true
    };
    config.fieldConfig.forcedE[0].f = ([](const std::array<T,2>& arr){return -arr[0];}); //Ex = -x
    config.fieldConfig.forcedE[1].f = ([](const std::array<T,2>& arr){return -arr[1];}); //Ey = -y
    config.fieldConfig.forcedB[2].f = ([](const std::array<T,2>& arr){return -0.01*arr[0]-0.01*arr[1];}); //Weak magnetic field
    return config;
}
