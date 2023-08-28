#ifndef CIE_UTILS_COMMAND_LINE_ABS_ARGUMENT_HPP
#define CIE_UTILS_COMMAND_LINE_ABS_ARGUMENT_HPP

// --- Utility Includes ---
#include "packages/commandline/inc/ArgParse.hpp"
#include "packages/macros/inc/typedefs.hpp"


namespace cie::utils::detail {


class AbsArgument
{
public:
    enum class Tag
    {
        Positional, // must be specified and has no default value

        Keyword, // identified by a flag/keyword, can be optional or required

        Flag, // optional binary value: it's either set or it isn't (has default)

        Invalid // invalid syntax
    }; // enum class Tag

    CIE_DEFINE_CLASS_POINTERS(AbsArgument)

public:
    AbsArgument(ArgParse::Key&& r_key);

    AbsArgument(const ArgParse::Key& r_key);

    CIE_DEFINE_CLASS_DEFAULTS(AbsArgument)

    virtual ~AbsArgument() = default;

    virtual bool matchesKey(const ArgParse::KeyView& r_key) const = 0;

    virtual bool matchesKey(const ArgParse::ValueView& r_key) const = 0;

    const ArgParse::Key& key() const;

private:
    ArgParse::Key _key;
};


} // namespace cie::utils::detail


#endif