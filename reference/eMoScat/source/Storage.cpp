#include <stdio.h>
#include <string>
#include <fstream>

#include "Storage.h"

namespace QSCAT
{
    bool BinaryStorageInterface::save_binary(const char* name) const
    {
        std::ofstream file;
        file.open(name, std::ios::out | std::ios::binary);
        bool state = save_bin_body(file);
        file.close();
        return state;
    }
    bool BinaryStorageInterface::save_binary(const std::string& name) const
    {
        std::ofstream file;
        file.open(name.c_str(), std::ios::out | std::ios::binary);
        bool state = save_bin_body(file);
        file.close();
        return state;
    }
    bool BinaryStorageInterface::save_binary(std::ofstream& file) const
    {
        return save_bin_body(file);
    }
    bool BinaryStorageInterface::read_binary(const char* name)
    {
        std::ifstream file;
        file.open(name, std::ios::in | std::ios::binary);
        bool state = read_bin_body(file);
        file.close();
        return state;
    }
    bool BinaryStorageInterface::read_binary(const std::string& name)
    {
        std::ifstream file;
        file.open(name.c_str(), std::ios::in | std::ios::binary);
        bool state = read_bin_body(file);
        file.close();
        return state;
    }
    bool BinaryStorageInterface::read_binary(std::ifstream& file)
    {
        return read_bin_body(file);
    }
}
