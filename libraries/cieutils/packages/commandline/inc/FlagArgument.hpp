#ifndef CIE_UTILS_COMMAND_LINE_FLAG_ARGUMENT_HPP
#define CIE_UTILS_COMMAND_LINE_FLAG_ARGUMENT_HPP

// --- Utiltiy Includes ---
#include "packages/commandline/inc/AbsArgument.hpp"

// --- STL Includes ---
#include <regex>


namespace cie::utils::detail {


class FlagArgument final : public AbsArgument
{
public:
    static const AbsArgument::Tag tag = AbsArgument::Tag::Flag;

public:
    FlagArgument(const ArgParse::Key& r_key);

    FlagArgument(FlagArgument&& r_rhs) = default;

    FlagArgument(const FlagArgument& r_rhs) = default;

    FlagArgument& operator=(FlagArgument&& r_rhs) = default;

    FlagArgument& operator=(const FlagArgument& r_rhs) = default;

    /// Match either the whole key (including leading the "--") or one character in a compound flag
    virtual bool matchesKey(const ArgParse::KeyView& r_key) const override;

    /// Match either the whole key (including leading the "--") or one character in a compound flag
    virtual bool matchesKey(const ArgParse::ValueView& r_key) const override;

private:
    std::regex _regex;
}; // class FlagArgument


} // namespace cie::utils::detail


#endif