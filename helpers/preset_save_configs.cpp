template<typename T, unsigned int Nd, unsigned int Nv>
typename Config<T,Nd,Nv>::SaveConfig preset_save_configs::Energies(T dt){
    typename Config<T,Nd,Nv>::SaveConfig saveConfig{false, false, false, false, true, false, false, true, false, false, true, dt, "outputs/"};
    return saveConfig;
}

template<typename T, unsigned int Nd, unsigned int Nv>
typename Config<T,Nd,Nv>::SaveConfig preset_save_configs::None(){
    typename Config<T,Nd,Nv>::SaveConfig saveConfig{false, false, false, false, false, false, false, false, false, false, false, 10e8, "outputs/"};
    return saveConfig;
}