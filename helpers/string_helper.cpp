template<typename T>
std::string string_helper::toString(T val){
    std::stringstream ss;
    try{
        ss << val;
        return ss.str();
    }catch(const std::exception& exception){
        throw std::invalid_argument("toString called with invalid input");
    }
}