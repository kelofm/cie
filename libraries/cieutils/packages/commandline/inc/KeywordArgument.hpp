#ifndef CIE_UTILS_COMMAND_LINE_KEYWORD_ARGUMENT_HPP
#define CIE_UTILS_COMMAND_LINE_KEYWORD_ARGUMENT_HPP

// --- Utility Includes ---
#include "packages/commandline/inc/AbsArgument.hpp"


namespace cie::utils::detail {


class KeywordArgument final : public AbsArgument
{
public:
    static const AbsArgument::Tag tag = AbsArgument::Tag::Keyword;

public:
    KeywordArgument(const ArgParse::Key& r_key);

    KeywordArgument(KeywordArgument&& r_rhs) = default;

    KeywordArgument(const KeywordArgument& r_rhs) = default;

    KeywordArgument& operator=(KeywordArgument&& r_rhs) = default;

    KeywordArgument& operator=(const KeywordArgument& r_rhs) = default;

    /// Matches only if the full input is identical to the stored key (including leading '-' characters)
    virtual bool matchesKey(const ArgParse::KeyView& r_key) const override;

    /// Matches only if the full input is identical to the stored key (including leading '-' characters)
    virtual bool matchesKey(const ArgParse::ValueView& r_key) const override;

private:
    ArgParse::Key _key;
}; // class KeywordArgument


} // namespace cie::utils::detail


#endif