#pragma once

// --- STL Includes ---
#include <concepts>


namespace cie::linalg {


template <class T>
concept LinalgSpaceLike
= requires (
    const T space,
    typename T::Vector vector,
    const typename T::Vector constVector,
    typename T::VectorView view,
    typename T::ConstVectorView constView) {
        typename T::Value;
        typename T::Vector;
        typename T::VectorView;
        typename T::ConstVectorView;

        {space.view(vector)}                                -> std::same_as<typename T::VectorView>;
        {space.view(constVector)}                           -> std::same_as<typename T::ConstVectorView>;
        {space.size(view)}                                  -> std::same_as<std::size_t>;
        {space.makeVector(std::size_t())}                   -> std::same_as<typename T::Vector>;

        {space.innerProduct(constView,constView)}           -> std::same_as<typename T::Value>;
        {space.scale(view, constView, typename T::Value())} -> std::same_as<void>;
        {space.scale(view, typename T::Value())}            -> std::same_as<void>;
        {space.add(view, constView, typename T::Value())}   -> std::same_as<void>;
        {space.assign(view, constView)}                     -> std::same_as<void>;
        {space.fill(view, typename T::Value())}             -> std::same_as<void>;
}; // concept LinalgSpaceLike


} // namespace cie::linalg
