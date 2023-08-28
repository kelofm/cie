// --- Utility Includes ---
#include "packages/commandline/inc/PositionalArgument.hpp"
#include "packages/exceptions/inc/exception.hpp"
#include "packages/macros/inc/exceptions.hpp"


namespace cie::utils::detail {


bool PositionalArgument::matchesKey(const ArgParse::KeyView&) const
{
    return false;
}


bool PositionalArgument::matchesKey(const ArgParse::ValueView&) const
{
    return false;
}


} // namespace cie::utils::detail