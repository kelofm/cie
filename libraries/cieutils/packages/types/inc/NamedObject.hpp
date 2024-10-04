#ifndef CIE_UTILS_NAMED_OBJECT_HPP
#define CIE_UTILS_NAMED_OBJECT_HPP

// --- STL Includes ---
#include <string>


namespace cie::utils {


/// @ingroup cieutils
class NamedObject
{
public:
    NamedObject(std::string&& r_name);

    NamedObject(const std::string& r_name);

    NamedObject(NamedObject&& r_rhs) = default;

    NamedObject(const NamedObject& r_rhs) = default;

    NamedObject& operator=(NamedObject&& r_rhs) = default;

    NamedObject& operator=(const NamedObject& r_rhs) = default;

    const std::string& name() const
    {return _name;}

    void setName(const std::string& r_name);

    virtual void setName(std::string::const_iterator begin,
                         std::string::const_iterator end);

    virtual ~NamedObject() {}

private:
    std::string _name;
};


} // namespace cie::utils

#endif