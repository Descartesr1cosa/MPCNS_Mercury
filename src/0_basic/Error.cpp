#include "0_basic/Error.h"
#include <cstdlib>

[[noreturn]] void ERROR::Abort(std::string Message)
{
    std::cout << "Error ! !\t" << Message << "\n"
              << std::flush;
    std::exit(-1);
}
