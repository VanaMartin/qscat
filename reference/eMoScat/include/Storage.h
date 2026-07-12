#ifndef INCLUDE_QSCAT_STORAGE_H_
#define INCLUDE_QSCAT_STORAGE_H_

#include <stdio.h>
#include <string>
#include <fstream>
namespace QSCAT
{
/** \addtogroup Interface
* @{ */

/// Simple binary I/O interface
/*!
    General interface to be inherited for all classes with storage
    capabilities. To successfully inherit the storage interface one must
    define protected "save_bin_body" and "read_bin_body" methods.
*/
class BinaryStorageInterface
{
 protected:
    /// Internal save to binary stream helper
    virtual bool save_bin_body(std::ofstream& file) const = 0;
    /// Internal read from bianry stream helper
    virtual bool read_bin_body(std::ifstream& file) = 0;
 public:
    /// saves class binary image into a separate file
    /*!
        To be used to store the instantion in file of a given name
        Arguments: name .. constant string of characters, given name with relative path
        Returns true on success false on fail
    */
    bool save_binary(const char *name) const;

    /// saves class binary image into a separate file
    /*!
        To be used to store the instantion in file of a given name
        Arguments: name .. constant string, given name with relative path
        Returns true on success false on fail
    */
    bool save_binary(const std::string& name) const;

    /// saves class binary image into binary stream
    /*!
        To be used for storing with other objects in a stream. Stream provider is responsible
        for the appropriate ordering within the stream
        Arguments: file .. std::ofstream, must be opened
        Returns true on success false on fail
    */
    bool save_binary(std::ofstream& file) const;

    /// read class binary image from a separate file
    /*!
        To be used to read a single instantion from a file
        Arguments: name .. const string of characters containing the file name
        Returns true on success false on fail
    */
    bool read_binary(const char *name);

    /// read class binary image from a separate file
    /*!
        To be used to read a single instantion from a file
        Arguments: name .. const string of characters containing the file name
        Returns true on success false on fail
    */
    bool read_binary(const std::string& name);

    /// read class binary image from binary stream
    /*!
        To be used for reading the instantion from stream of binary data. Sream provider is responsible
        for the appropriate ordering within the stream.
        Arguments: file .. std::ifstream, must be opened
        Returns true on success false on fail
    */
    bool read_binary(std::ifstream& file);
};

/** @} */
} // namespace QSCAT

#endif // INCLUDE_QSCAT_STORAGE_H_
