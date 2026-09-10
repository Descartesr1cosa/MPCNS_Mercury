#include "MarsDipoleCache.h"
#include <fstream>
#include <filesystem>
#include <cmath>
#include <cstring>

namespace MARS
{
namespace
{
constexpr char magic[8]={'M','A','R','S','B','A','D','1'};
std::uint64_t checksum(const std::vector<double> &v)
{DipoleFingerprint h;h.Bytes(v.data(),v.size()*sizeof(double));return h.value;}
}
bool ReadDipoleCache(const std::string &path,std::uint64_t signature,std::vector<double> &values)
{
    std::ifstream in(path,std::ios::binary);
    char tag[8]{};std::uint64_t saved=0,count=0,hash=0;
    in.read(tag,8);in.read(reinterpret_cast<char*>(&saved),8);
    in.read(reinterpret_cast<char*>(&count),8);in.read(reinterpret_cast<char*>(&hash),8);
    if(!in || std::memcmp(tag,magic,8)!=0 || saved!=signature || count!=values.size())return false;
    std::vector<double> candidate(values.size());
    in.read(reinterpret_cast<char*>(candidate.data()),candidate.size()*sizeof(double));
    if(!in || in.peek()!=std::char_traits<char>::eof() || checksum(candidate)!=hash)return false;
    for(double v:candidate)if(!std::isfinite(v))return false;
    values.swap(candidate);return true;
}
bool WriteDipoleCache(const std::string &path,std::uint64_t signature,const std::vector<double> &values)
{
    try
    {
        std::filesystem::create_directories(std::filesystem::path(path).parent_path());
        const std::string temporary=path+".tmp";
        std::ofstream out(temporary,std::ios::binary|std::ios::trunc);
        const std::uint64_t count=values.size(),hash=checksum(values);
        out.write(magic,8);out.write(reinterpret_cast<const char*>(&signature),8);
        out.write(reinterpret_cast<const char*>(&count),8);out.write(reinterpret_cast<const char*>(&hash),8);
        out.write(reinterpret_cast<const char*>(values.data()),values.size()*sizeof(double));
        out.close();if(!out)return false;
        std::filesystem::rename(temporary,path);return true;
    }
    catch(const std::filesystem::filesystem_error &){return false;}
}
}
