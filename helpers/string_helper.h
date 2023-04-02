#ifndef PIC_SEMI_IMPLICIT_STRING_HELPER_H
#define PIC_SEMI_IMPLICIT_STRING_HELPER_H

#include <string>
#include <sstream>

namespace string_helper {
    template<typename T>
    std::string toString(T val);
}

#include "string_helper.cpp"

#endif //PIC_SEMI_IMPLICIT_STRING_HELPER_H
