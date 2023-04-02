template <typename T>
typename Config<T,1,1>::SpeciesConfig preset_species::Uniform1D1V(unsigned int Np, T m, T q, T Lx, unsigned int Nx, T Kb, T T0){
    return typename Config<T,1,1>::SpeciesConfig{
        Np,
        m,
        q,
        preset_distributions::Uniform<double,1>(1),
        Grid<T,1>{std::array<typename Grid<T,1>::Dim,1>{{{0,Lx,Nx}}}},
        preset_distributions::Boltzmann<double,1>(m,Kb,T0),
        Grid<T,1>{std::array<typename Grid<T,1>::Dim,1>{{
            {-math_helper::boltzmannBounds(m,Kb,T0),math_helper::boltzmannBounds(m,Kb,T0),1000}
        }}},
        false,
        "",
        false,
        ""
    };
}

template <typename T>
typename Config<T,2,3>::SpeciesConfig preset_species::Uniform2D3V(unsigned int Np, T m, T q, T Lx, T Ly, unsigned int Nx, unsigned int Ny, T Kb, T T0){
    return typename Config<T,2,3>::SpeciesConfig{
            Np,
            m,
            q,
            preset_distributions::Uniform<double,2>(1),
            Grid<T,2>{std::array<typename Grid<T,2>::Dim,2>{{{0,Lx,Nx},{0,Ly,Ny}}}},
            preset_distributions::Boltzmann<double,3>(m,Kb,T0),
            Grid<T,3>{std::array<typename Grid<T,3>::Dim,3>{{
                {-math_helper::boltzmannBounds(m,Kb,T0),math_helper::boltzmannBounds(m,Kb,T0),100},
                {-math_helper::boltzmannBounds(m,Kb,T0),math_helper::boltzmannBounds(m,Kb,T0),100},
                {-math_helper::boltzmannBounds(m,Kb,T0),math_helper::boltzmannBounds(m,Kb,T0),100}
            }}},
            false,
            "",
            false,
            ""
    };
}
