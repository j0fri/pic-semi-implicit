template <typename T>
typename Config<T,1,1>::FieldConfig preset_fields::Default1D1V(T Lx, unsigned int Nx, T c, T e0) {
    return typename Config<T,1,1>::FieldConfig{
        Grid<T,1>{std::array<typename Grid<T,1>::Dim,1>{{{0,Lx,Nx}}}},
        c,
        e0,
        {preset_distributions::Uniform<T,1>(0)},
        {preset_distributions::Uniform<T,1>(0)},
        false,
        false,
        true
    };
}


template <typename T>
typename Config<T,2,3>::FieldConfig preset_fields::Default2D3V(T Lx, T Ly, unsigned int Nx, unsigned int Ny, T c, T e0) {
    return typename Config<T,2,3>::FieldConfig{
            Grid<T,2>{std::array<typename Grid<T,2>::Dim,2>{{{0,Lx,Nx},{0,Ly,Ny}}}},
            c,
            e0,
            {preset_distributions::Uniform<T,2>(0),preset_distributions::Uniform<T,2>(0),preset_distributions::Uniform<T,2>(0)},
            {preset_distributions::Uniform<T,2>(0),preset_distributions::Uniform<T,2>(0),preset_distributions::Uniform<T,2>(0)},
            false,
            false,
            true
    };
}


