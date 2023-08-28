#ifndef CIE_UTILS_PARAMETER_PACK_ARGUMENT_MAPPER_IMPL_HPP
#define CIE_UTILS_PARAMETER_PACK_ARGUMENT_MAPPER_IMPL_HPP

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"


namespace cie::ct {


namespace detail {


template <Size I, class ...TDefaults, class ...TSpecified>
requires (0 == I)
void insertSpecified(std::tuple<TDefaults...>& r_defaults, const std::tuple<TSpecified...>& r_specified)
{
    CIE_BEGIN_EXCEPTION_TRACING

    using Current = typename std::decay<typename std::tuple_element<I,std::tuple<TSpecified...>>::type>::type;
    std::get<Current>(r_defaults) = std::get<I>(r_specified);

    CIE_END_EXCEPTION_TRACING
}


template <Size I, class ...TDefaults, class ...TSpecified>
requires (0 < I)
void insertSpecified(std::tuple<TDefaults...>& r_defaults, const std::tuple<TSpecified...>& r_specified)
{
    CIE_BEGIN_EXCEPTION_TRACING

    using Current = typename std::decay<typename std::tuple_element<I,std::tuple<TSpecified...>>::type>::type;
    std::get<Current>(r_defaults) = std::get<I>(r_specified);
    insertSpecified<I-1>(r_defaults, r_specified);

    CIE_END_EXCEPTION_TRACING
}


} // namespace detail


template <class ...TDefaults>
ArgumentMapper<TDefaults...>::ArgumentMapper(TDefaults... r_defaults)
    : _defaults(std::forward<TDefaults>(r_defaults)...)
{
}


template <class ...TDefaults>
ArgumentMapper<TDefaults...>::ArgumentMapper(Tuple&& r_defaults)
    : _defaults(std::move(r_defaults))
{
}


template <class ...TDefaults>
template <class TReturn, class TFunction, class ...TSpecified>
requires (!Empty<TSpecified...>)
inline TReturn
ArgumentMapper<TDefaults...>::mapImpl(TFunction&& r_function, std::tuple<TSpecified...>&& r_specified)
{
    Tuple passed = _defaults;
    detail::insertSpecified<ct::PackSize<TSpecified...>-1>(passed, std::move(r_specified));
    return std::apply(std::forward<TFunction>(r_function), std::move(passed));
}


template <class ...TDefaults>
template <class TReturn, class TFunction, class ...TSpecified>
requires Empty<TSpecified...>
inline TReturn
ArgumentMapper<TDefaults...>::mapImpl(TFunction&& r_function, std::tuple<TSpecified...>&&)
{
    return std::apply(std::forward<TFunction>(r_function), Tuple(_defaults));
}


template <class ...TDefaults>
template <class TReturn, class TFunction, class ...TSpecified>
requires concepts::AppliableWith<TFunction,TDefaults...>
inline TReturn
ArgumentMapper<TDefaults...>::map(TFunction&& r_function, TSpecified&&... r_specified)
{
    CIE_BEGIN_EXCEPTION_TRACING

    return this->mapImpl<TReturn>(
        std::forward<TFunction>(r_function),
        std::tuple<TSpecified&&...>(std::forward<TSpecified>(r_specified)...)
    );

    CIE_END_EXCEPTION_TRACING
}


template <class ...TDefaults>
template <class TReturn, class TFunction, class TTuple>
inline TReturn
ArgumentMapper<TDefaults...>::mapTuple(TFunction&& r_function, TTuple&& r_specified)
{
    CIE_BEGIN_EXCEPTION_TRACING

    return this->mapImpl<TReturn>(
        std::forward<TFunction>(r_function),
        std::move(r_specified)
    );

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::ct


#endif