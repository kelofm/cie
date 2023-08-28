// --- Utility Includes ---
#include "packages/commandline/inc/KeywordArgument.hpp"
#include "packages/commandline/inc/KeywordParser.hpp"
#include "packages/stl_extension/inc/zip.hpp"
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/checks.hpp"




namespace cie::utils::detail {


KeywordArgument::KeywordArgument(const ArgParse::Key& r_key)
    : AbsArgument(r_key),
      _key(KeywordParser::parse(r_key))
{
}


bool KeywordArgument::matchesKey(const ArgParse::KeyView& r_key) const
{
    return r_key == _key;
}


bool KeywordArgument::matchesKey(const ArgParse::ValueView& r_key) const
{
    if (r_key.size() != _key.size())
        return false;

    // Cannot zip because ranges have no value_type >.<'
    //for (auto [value, reference] : zip(r_key, _key))
    //    if (value != reference)
    //        return false;

    auto it = r_key.begin();
    const auto it_end = r_key.end();
    auto it_ref = _key.begin();
    for (; it!=it_end; ++it, ++it_ref)
        if ((*it) != (*it_ref))
            return false;

    return true;
}


} // namespace cie::utils::detail
