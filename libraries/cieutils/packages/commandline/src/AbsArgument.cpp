// --- Utility Includes ---
#include "packages/commandline/inc/AbsArgument.hpp"


namespace cie::utils::detail {


AbsArgument::AbsArgument(ArgParse::Key&& r_key)
    : _key(std::move(r_key))
{
}


AbsArgument::AbsArgument(const ArgParse::Key& r_key)
    : _key(r_key)
{
}


const ArgParse::Key& AbsArgument::key() const
{
    return _key;
}


} // namespace cie::utils::detail
