template <typename T>
typename Config<T,1,1>::SpeciesConfig preset_species::Uniform1D1V(unsigned int Np, T m, T q, T Lx, unsigned int Nx, T Kb, T T0){
    return typename Config<T,1,1>::SpeciesConfig{
        Np,
        m,
        q,
        preset_distributions::Uniform<double,1>(1),
        Grid<T,1>{std::array<typename Grid<T,1>::Dim,1>{{{0,Lx,Nx}}}},
        preset_distributions::Boltzmann<double,1>(1,Kb,T0),
        Grid<T,1>{std::array<typename Grid<T,1>::Dim,1>{{{0,Lx,1000}}}},
        false,
        "",
        false,
        ""
    };
}
