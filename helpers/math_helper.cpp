template <typename T>
std::vector<std::vector<T>> math_helper::cartesian(const std::vector<std::vector<T>>& s1, const std::vector<T>& s2){
	int height1 = s1.size();
	int height2 = s2.size();
	int width = s1[0].size();
	std::vector<std::vector<T>> out(height1*height2, std::vector<T>(width+1,0));
	int index = 0;
	for(int i = 0; i < height1; ++i){
		for(int j = 0; j < height2; ++j){
			std::copy(s1[i].begin(), s1[i].end(), out[index].begin());
			out[index].back() = s2[j];
			++index;
		}
	} 
	return out;
}

template <typename T>
std::vector<std::vector<T>> math_helper::repeatedCartesian(const std::vector<std::vector<T>>& sets){
	int n = sets.size();
	const std::vector<T>& first = sets[0];
	std::vector<std::vector<T>> out(first.size(),std::vector<T>(1));
	for(int i = 0; i < (int)first.size(); ++i){
		out[i][0] = first[i];
	}
	for(int i = 1; i < n; ++i){
		out = math_helper::cartesian(out, sets[i]);
	}
	return out;
}

template <typename T>
std::vector<T> math_helper::linspace(T x1, T x2, unsigned int n){
	std::vector<T> out(n);
	T dx = (x2-x1)/(n-1);
	out.front() = x1;
	for(unsigned int i = 1; i < n-1; ++i){
		out[i] = x1 + i * dx;
	}
	out.back() = x2;
	return out;
}

//TODO: consider matrix multiplication as addition of columns, might improve cache locality
template<typename T>
void math_helper::gemv(unsigned int M, unsigned int N, T alpha, const T* A, unsigned int lda, const T* x,
                       unsigned int incx, T* y, unsigned int incy){
    T temp;
    for(unsigned int i = 0; i < M; ++i){
        temp = 0;
        for(unsigned int j = 0; j < N; ++j){
            temp += A[i+lda*j] * x[j*incx];
        }
        y[i*incy] = alpha*temp;
    }
}

template <typename T>
T math_helper::boltzmannBounds(T m, T Kb, T T0) {
    return (T)4.2 * std::sqrt(Kb*T0/m);
}
