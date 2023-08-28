#ifndef CIE_UTILS_MAKE_SHARED_FROM_TUPLE_IMPL_HPP
#define CIE_UTILS_MAKE_SHARED_FROM_TUPLE_IMPL_HPP

// --- Internal Includes ---
#include "packages/macros/inc/exceptions.hpp"


namespace cie::utils {


namespace detail {

template <class T, class ...Args>
struct make_shared_from_tuple_helper
{
    make_shared_from_tuple_helper(std::tuple<Args...>& args) :
        _args(args) {}

    std::shared_ptr<T> make()
    {
        return std::apply(std::make_shared<T,std::add_lvalue_reference_t<Args>...>, _args);
    }

private:
    std::tuple<Args...>& _args;
};

} // namespace detail


template <class T, class ...Args>
std::shared_ptr<T> make_shared_from_tuple(std::tuple<Args...>& args)
{
    CIE_BEGIN_EXCEPTION_TRACING
    return detail::make_shared_from_tuple_helper<T,Args...>(args).make();
    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::utils


#endif