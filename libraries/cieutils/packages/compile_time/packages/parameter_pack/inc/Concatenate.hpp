#ifndef CIE_UTILS_PARAMETER_PACK_CONCATENATE_HPP
#define CIE_UTILS_PARAMETER_PACK_CONCATENATE_HPP

// --- STL Includes ---
#include <tuple>


namespace cie::ct {


template <class ...Ts>
struct Concatenate
{
    template <class ...TTs>
    using Type = std::tuple<Ts...,TTs...>;
}; // struct Concatenate


template <class T>
struct ConcatenateTuple {};


template <class ...Ts>
struct ConcatenateTuple<std::tuple<Ts...>>
{
private:
    template <class T>
    struct TypeImpl {};

    template <class ...TTs>
    struct TypeImpl<std::tuple<TTs...>>
    {using Type = std::tuple<Ts...,TTs...>;};

public:
    template <class T>
    using Type = typename TypeImpl<T>::Type;
}; // struct ConcatenateTuple


} // namespace cie::ct


#endif
