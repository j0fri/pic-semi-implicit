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
Config<T,2,3> preset_configs::landau2D3VX(unsigned int Nx, unsigned int Ny){
    Config<T,2,3> config{
        std::vector<typename Config<T,2,3>::SpeciesConfig>{{
               preset_species::Uniform2D3V<T>(50000,1,-1,4*M_PI,1,Nx,Ny,1,1),
               preset_species::Uniform2D3V<T>(50000,2000,1,4*M_PI,1,Nx,Ny,1,1),
        }},
        preset_fields::Default2D3V<T>(4*M_PI,1,Nx,Ny,1,1),
        {10,0.1},
        {false, true, false, false, false, true, true, false, true, 1, "outputs/"},
        {{true,true},{}}, //Only periodic boundary conditions
        true,
        true
    };
    //Add perturbation in electrons
    config.speciesConfig[0].xDist += preset_distributions::Sin<T,2>(0,0.2,(T)1.0/2,0);

    return config;
}

template <typename T>
Config<T,2,3> preset_configs::constAccelerationX(){
    Config<T,2,3> config{
            std::vector<typename Config<T,2,3>::SpeciesConfig>{{
                   preset_species::TopHat2D3VX<double>(50000,1,1,-0.2,0.2,-0.2,0.2)
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
                   preset_species::TopHat2D3VX<double>(50000,1,1,-0.2,0.2,0.3,0.7)
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
