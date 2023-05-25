template<typename T, unsigned int Nd, unsigned int Nv>
typename Config<T,Nd,Nv>::SaveConfig preset_save_configs::Energies(T dt){
    typename Config<T,Nd,Nv>::SaveConfig saveConfig{false, false, false, false, true, false, false, true, false, false, dt, "outputs/"};
    return saveConfig;
}