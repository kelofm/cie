#ifndef CIE_UTILS_RANGES_TRANSFORM_VIEW_HPP
#define CIE_UTILS_RANGES_TRANSFORM_VIEW_HPP

// --- STL Includes ---
#include <ranges>
#include <functional>


namespace cie::utils {


namespace detail {

template <class ContainerType, class OutputType>
struct TransformViewTraits
{
    using stored_type         = typename ContainerType::value_type;
    using value_type          = OutputType;
    using container_type      = ContainerType;
    using container_view_type = std::ranges::views::all_t<const container_type&>;
    using transform_type      = std::function<OutputType(const typename ContainerType::value_type&)>;

    using size_type           = typename ContainerType::size_type;
    using iterator            = typename ContainerType::iterator;
    using const_iterator      = typename ContainerType::const_iterator;
};

} // namespace detail


/// @addtogroup cieutils
template <class ContainerType, class OutputType>
class TransformView :
    public std::ranges::transform_view< typename detail::TransformViewTraits<ContainerType,OutputType>::container_view_type,
                                        typename detail::TransformViewTraits<ContainerType,OutputType>::transform_type >,
    public detail::TransformViewTraits<ContainerType,OutputType>
{
protected:
    using view_type = std::ranges::transform_view< typename detail::TransformViewTraits<ContainerType,OutputType>::container_view_type,
                                                   typename detail::TransformViewTraits<ContainerType,OutputType>::transform_type >;

public:
    template <class ...Args>
    TransformView( Args&&... args ) :
        view_type( std::forward<Args>(args)... )
    {}

    TransformView& operator=(TransformView&& r_rhs) noexcept = default;

    TransformView& operator=(const TransformView& r_rhs) = default;
};


/// @addtogroup cieutils
template <class OutputType, class ContainerType>
TransformView<ContainerType,OutputType> makeTransformView( const ContainerType& r_container,
                                                           typename TransformView<ContainerType,OutputType>::transform_type transformFunctor )
{
    return TransformView<ContainerType,OutputType>( r_container, transformFunctor );
}


} // namespace cie::utils


#endif