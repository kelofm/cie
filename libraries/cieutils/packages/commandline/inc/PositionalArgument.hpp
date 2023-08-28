#ifndef CIE_UTILS_COMMAND_LINE_POSITIONAL_ARGUMENT_HPP
#define CIE_UTILS_COMMAND_LINE_POSITIONAL_ARGUMENT_HPP

// --- Utility Includes ---
#include "packages/commandline/inc/AbsArgument.hpp"


namespace cie::utils::detail {


class PositionalArgument final : public AbsArgument
{
public:
    static const AbsArgument::Tag tag = AbsArgument::Tag::Positional;

public:
    using AbsArgument::AbsArgument;

    /// Never matches
    virtual bool matchesKey(const ArgParse::KeyView& r_key) const override;

    /// Never matches
    virtual bool matchesKey(const ArgParse::ValueView& r_key) const override;
}; // class PositionalArgument


} // namespace cie::utils::detail


#endif