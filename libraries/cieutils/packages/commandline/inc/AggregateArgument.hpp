#ifndef CIE_UTILS_COMMAND_LINE_AGGREGATE_ARGUMENT_HPP
#define CIE_UTILS_COMMAND_LINE_AGGREGATE_ARGUMENT_HPP

// --- Utility Includes ---
#include "packages/commandline/inc/AbsArgument.hpp"
#include "packages/macros/inc/typedefs.hpp"
#include "packages/types/inc/IDObject.hpp"
#include "packages/types/inc/NamedObject.hpp"


namespace cie::utils::detail {


class AbsAggregateArgument : public AbsArgument
{
public:
    CIE_DEFINE_CLASS_POINTERS(AbsAggregateArgument)

    using Map = ArgParse::Results::Map;

public:
    virtual void parseValues(const ArgParse::ValueViewContainer& r_input, Map& r_outputMap) const = 0;

    virtual const ArgParse::DefaultValue& defaultValue() const = 0;

    virtual bool isOptional() const = 0;

    virtual ArgParse::ArgumentCount nArgs() const = 0;

    virtual const std::string& docString() const = 0;

    virtual AbsArgument::Tag tag() const = 0;

    virtual const std::string& name() const = 0;

    virtual ArgParse::KeyContainer keys() const = 0;

    static Tag parseTag(const ArgParse::Key& r_key, ArgParse::ArgumentCount nArgs);
}; // class AbsAggregateArgument


template <class TArgument>
class AggregateArgument final : public AbsAggregateArgument, public IDObject<Size>, public NamedObject
{
public:
    using ArgumentContainer = DynamicArray<TArgument>;

public:
    AggregateArgument(ArgumentContainer&& r_arguments,
                      Size id,
                      bool isOptional,
                      ArgParse::DefaultValue&& r_defaultValue,
                      ArgParse::ArgumentCount nArgs,
                      const ArgParse::Validator& r_validator,
                      std::string&& r_docString,
                      std::string&& r_displayName);

    AggregateArgument(AggregateArgument&& r_rhs) = default;

    AggregateArgument(const AggregateArgument& r_rhs) = default;

    AggregateArgument& operator=(AggregateArgument&& r_rhs) = default;

    AggregateArgument& operator=(const AggregateArgument& r_rhs) = default;

    template <class ...TArgs>
    void add(TArgs&&... r_args)
    {_arguments.emplace_back(std::forward<TArgs>(r_args)...);}

    /// Return true if any of the stored arguments matches the provided key
    virtual bool matchesKey(const ArgParse::KeyView& r_key) const override;

    /// Return true if any of the stored arguments matches the provided key
    virtual bool matchesKey(const ArgParse::ValueView& r_key) const override;

    virtual void parseValues(const ArgParse::ValueViewContainer& r_input, Map& r_outputMap) const override;

    virtual const ArgParse::DefaultValue& defaultValue() const override;

    virtual bool isOptional() const override;

    virtual ArgParse::ArgumentCount nArgs() const override;

    virtual const std::string& docString() const override;

    virtual AbsArgument::Tag tag() const override;

    /// Shadows NamedObject::name
    virtual const std::string& name() const override;

    virtual ArgParse::KeyContainer keys() const override;

private:
    bool _isOptional;

    ArgParse::DefaultValue _defaultValue;

    ArgParse::ArgumentCount _nArgs;

    ArgParse::Validator _validator;

    std::string _docString;

    ArgumentContainer _arguments;
}; // class AggregateArgument


OutputStream& operator<<(OutputStream& r_stream, const AbsAggregateArgument& r_argument);


std::ostream& operator<<(std::ostream& r_stream, const AbsAggregateArgument& r_argument);


} // namespace cie::utils::detail


#endif