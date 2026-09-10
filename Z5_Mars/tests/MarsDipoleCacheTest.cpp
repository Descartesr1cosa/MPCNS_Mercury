#include "MarsDipoleCache.h"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <limits>

void check(bool value){if(!value)throw std::runtime_error("dipole cache check failed");}
int main()
{
    const auto dir=std::filesystem::temp_directory_path()/
        ("mars-cache-test-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    const auto path=(dir/"dipole.bin").string();
    const std::vector<double> source{1,-2,3,0};std::vector<double> data(4,99);
    check(!MARS::ReadDipoleCache(path,12,data));
    check(MARS::WriteDipoleCache(path,12,source));
    check(!MARS::ReadDipoleCache(path,13,data));check(data[0]==99);
    std::vector<double> wrong_size(3);check(!MARS::ReadDipoleCache(path,12,wrong_size));
    check(MARS::ReadDipoleCache(path,12,data));check(data==source);
    {std::ofstream out(path,std::ios::binary|std::ios::app);out.put('x');}
    check(!MARS::ReadDipoleCache(path,12,data));
    check(MARS::WriteDipoleCache(path,12,source));
    {std::fstream out(path,std::ios::binary|std::ios::in|std::ios::out);out.seekp(32);out.put('x');}
    check(!MARS::ReadDipoleCache(path,12,data));check(data==source);
    std::filesystem::resize_file(path,12);check(!MARS::ReadDipoleCache(path,12,data));
    check(MARS::WriteDipoleCache(path,12,{1,2,3,std::numeric_limits<double>::infinity()}));
    check(!MARS::ReadDipoleCache(path,12,data));
    std::filesystem::remove_all(dir);
}
