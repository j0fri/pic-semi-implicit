template <typename T>
typename Config<T,1,1>::SpeciesConfig preset_species::Uniform1D1V(unsigned int Np, T m, T q, T Lx, unsigned int Nx, T Kb, T T0){
    return typename Config<T,1,1>::SpeciesConfig{
        Np,
        m,
        q,
        preset_distributions::Uniform<T,1>(1),
        Grid<T,1>{std::array<typename Grid<T,1>::Dim,1>{{{0,Lx,Nx}}}},
        preset_distributions::Boltzmann<T,1>(m,Kb,T0),
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
            preset_distributions::Uniform<T,2>(1),
            Grid<T,2>{std::array<typename Grid<T,2>::Dim,2>{{{0,Lx,Nx},{0,Ly,Ny}}}},
            preset_distributions::Boltzmann<T,3>(m,Kb,T0),
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

template <typename T>
typename Config<T,2,3>::SpeciesConfig preset_species::Uniform2D3VCustomBounds(unsigned int Np, T m, T q, T Lx, T Ly, unsigned int Nx, unsigned int Ny, T Kb, T T0, T zScore){
    return typename Config<T,2,3>::SpeciesConfig{
            Np,
            m,
            q,
            preset_distributions::Uniform<T,2>(1),
            Grid<T,2>{std::array<typename Grid<T,2>::Dim,2>{{{0,Lx,Nx},{0,Ly,Ny}}}},
            preset_distributions::Boltzmann<T,3>(m,Kb,T0),
            Grid<T,3>{std::array<typename Grid<T,3>::Dim,3>{{
                    {-math_helper::boltzmannBounds(m,Kb,T0,zScore),math_helper::boltzmannBounds(m,Kb,T0,zScore),100},
                    {-math_helper::boltzmannBounds(m,Kb,T0,zScore),math_helper::boltzmannBounds(m,Kb,T0,zScore),100},
                    {-math_helper::boltzmannBounds(m,Kb,T0,zScore),math_helper::boltzmannBounds(m,Kb,T0,zScore),100}
            }}},
            false,
            "",
            false,
            ""
    };
}

template <typename T>
typename Config<T,2,3>::SpeciesConfig preset_species::TopHat2D3V(unsigned int Np, T m, T q, T x1, T x2, T y1, T y2){
    return typename Config<T,2,3>::SpeciesConfig{
            Np,
            m,
            q,
            preset_distributions::Uniform<T,2>(1),
            Grid<T,2>{std::array<typename Grid<T,2>::Dim,2>{{{x1,x2,1},{y1,y2,1}}}},
            preset_distributions::Uniform<T,3>(1),
            Grid<T,3>{std::array<typename Grid<T,3>::Dim,3>{{
                    {0,0,1},{0,0,1},{0,0,1}
            }}}, //All velocities are zero
            false,
            "",
            false,
            ""
    };
}

template <typename T>
typename Config<T,2,3>::SpeciesConfig preset_species::TopHat2D3VBoltzmann(unsigned int Np, T m, T q, T x1, T x2, T y1, T y2, T Kb, T T0){
    return typename Config<T,2,3>::SpeciesConfig{
            Np,
            m,
            q,
            preset_distributions::Uniform<T,2>(1),
            Grid<T,2>{std::array<typename Grid<T,2>::Dim,2>{{{x1,x2,1},{y1,y2,1}}}},
            preset_distributions::Boltzmann<T,3>(m,Kb,T0),
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

template <typename T>
typename Config<T,2,3>::SpeciesConfig preset_species::fromFile(unsigned int Np, T m, T q, std::string posFileName, std::string velFileName){
    return typename Config<T,2,3>::SpeciesConfig{
            Np,
            m,
            q,
            preset_distributions::Uniform<T,2>(1),//No effect
            Grid<T,2>{std::array<typename Grid<T,2>::Dim,2>{{{0,0,1},{0,0,1}}}},//No effect
            preset_distributions::Uniform<T,3>(1),//No effect
            Grid<T,3>{std::array<typename Grid<T,3>::Dim,3>{{
                    {0,0,1},{0,0,1},{0,0,1}
            }}}, //No effect
            true,
            posFileName,
            true,
            velFileName
    };
}


