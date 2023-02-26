template <typename T>
Config<T,1,1> preset_configs::landau1D1V(){
     Config<double,1,1> config{
        std::vector<Config<double,1,1>::SpeciesConfig>{{
               preset_species::Uniform1D1V<double>(1000,1,-1,4*M_PI,30,1,1),
               preset_species::Uniform1D1V<double>(1000,1000,1,4*M_PI,30,1,1),
        }},
        preset_fields::Default1D1V<T>(4*M_PI,30,1,1),
        {10,0.01},
        {false, false, false, false, false, true, false, false, 0.5, "", "outputs/fieldOutput.txt"}
    };
    //Add perturbation in electrons
    config.speciesConfig[0].xDist += preset_distributions::Sin<T,1>(0,0.2,(T)1.0/2,0);
    return config;
}
