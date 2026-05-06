#pragma once
#include <iostream>
#include <string>

namespace ERROR
{
    [[noreturn]] void Abort(std::string Message);
} // namespace ERROR
