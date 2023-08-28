// --- Internal Includes ---
#include "packages/types/inc/NamedObject.hpp"


namespace cie::utils {


NamedObject::NamedObject(std::string&& r_name)
    : _name(std::move(r_name))
{
}


NamedObject::NamedObject(const std::string& r_name) :
    _name(r_name)
{
}


void NamedObject::setName(const std::string& r_name)
{
    this->setName(r_name.begin(), r_name.end());
}


void NamedObject::setName(std::string::const_iterator begin,
                          std::string::const_iterator end)
{
    _name = std::string(begin, end);
}


} // namespace cie::utils