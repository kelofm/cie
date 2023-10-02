#ifndef CIE_UTILS_PARAMETER_PACK_ARGUMENT_MAPPER_HPP
#define CIE_UTILS_PARAMETER_PACK_ARGUMENT_MAPPER_HPP

// --- Utility Includes ---
#include "packages/macros/inc/typedefs.hpp"
#include "packages/macros/inc/exceptions.hpp"
#include "packages/compile_time/packages/parameter_pack/inc/Size.hpp"
#include "packages/compile_time/packages/concepts/inc/functional.hpp"

// --- STL Includes ---
#include <tuple>


namespace cie::ct {


///@addtogroup cieutils
///@{

/** Utility for creating flexible interfaces for functions with many arguments.
 *
 *  @param TDefaults: default argument types
 *  @note @p TDefaults must be a unique set of types.
 */
template <class ...TDefaults>
class ArgumentMapper
{
private:
    using Tuple = std::tuple<TDefaults...>;

public:
    CIE_DEFINE_CLASS_DEFAULTS(ArgumentMapper)

    ArgumentMapper(TDefaults... r_defaults);

    ArgumentMapper(Tuple&& r_defaults);

    /// Call the provided function with the specified arguments, filling in the missing arguments with defaults.
    template <class TReturn, class TFunction, class ...TSpecified>
    requires concepts::AppliableWith<TFunction,TDefaults...>
    TReturn map(TFunction&& r_function, TSpecified&&... r_specified)
    {
        CIE_BEGIN_EXCEPTION_TRACING

        return this->mapImpl<TReturn>(
            std::forward<TFunction>(r_function),
            std::tuple<TSpecified&&...>(std::forward<TSpecified>(r_specified)...)
        );

        CIE_END_EXCEPTION_TRACING
    }

    template <class TReturn, class TFunction, class TTuple>
    TReturn mapTuple(TFunction&& r_function, TTuple&& r_specified);

private:
    template <class TReturn, class TFunction, class ...TSpecified>
    requires (!Empty<TSpecified...>)
    TReturn mapImpl(TFunction&& r_function, std::tuple<TSpecified...>&& r_specified);

    template <class TReturn, class TFunction, class ...TSpecified>
    requires Empty<TSpecified...>
    TReturn mapImpl(TFunction&& r_function, std::tuple<TSpecified...>&& r_specified);

private:
    Tuple _defaults;
};

///@}


} // namespace cie::ct

#include "packages/compile_time/packages/parameter_pack/impl/ArgumentMapper_impl.hpp"

#endif