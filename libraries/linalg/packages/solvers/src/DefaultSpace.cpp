// --- Linalg Includes ---
#include "packages/solvers/inc/DefaultSpace.hpp"

// --- Utility Includes ---
#include "packages/concurrency/inc/ParallelFor.hpp"
#include "packages/concurrency/inc/ThreadPoolBase.hpp"

// --- STL Includes ---
#include <numeric>
#include <algorithm>


namespace cie::linalg {


template <class T, TagLike TTag>
requires (std::is_same_v<TTag,tags::Serial> || std::is_same_v<TTag,tags::SMP>)
DefaultSpace<T,TTag>::DefaultSpace(Ref<mp::ThreadPoolBase> rThreads)
requires std::is_same_v<TTag,tags::SMP>
    : _maybeThreads(rThreads)
{}


template <class T, TagLike TTag>
requires (std::is_same_v<TTag,tags::Serial> || std::is_same_v<TTag,tags::SMP>)
typename DefaultSpace<T,TTag>::VectorView
DefaultSpace<T,TTag>::view(Ref<Vector> rVector) noexcept {
    return rVector;
}


template <class T, TagLike TTag>
requires (std::is_same_v<TTag,tags::Serial> || std::is_same_v<TTag,tags::SMP>)
typename DefaultSpace<T,TTag>::ConstVectorView
DefaultSpace<T,TTag>::view(Ref<const Vector> rVector) noexcept {
    return rVector;
}


template <class T, TagLike TTag>
requires (std::is_same_v<TTag,tags::Serial> || std::is_same_v<TTag,tags::SMP>)
std::size_t
DefaultSpace<T,TTag>::size(ConstVectorView view) noexcept {
    return view.size();
}


template <class T, TagLike TTag>
requires (std::is_same_v<TTag,tags::Serial> || std::is_same_v<TTag,tags::SMP>)
typename DefaultSpace<T,TTag>::Vector
DefaultSpace<T,TTag>::makeVector(std::size_t size) {
    return Vector(size);
}


template <class T, TagLike TTag>
requires (std::is_same_v<TTag,tags::Serial> || std::is_same_v<TTag,tags::SMP>)
typename DefaultSpace<T,TTag>::Value
DefaultSpace<T,TTag>::innerProduct(ConstVectorView left, ConstVectorView right) const {
    if (_maybeThreads.has_value()) {
        // Do a 2-stage reduction.
        auto loop = mp::ParallelFor<std::size_t>(_maybeThreads.value())
            .firstPrivate(static_cast<Value>(0));

        // Do the first stage in parallel.
        loop.execute(
            this->size(left),
            [left, right] (std::size_t iComponent, Ref<Value> rSum) -> void {
                rSum += left[iComponent] * right[iComponent];
            });

        // Do the second stage in serial.
        Value output = static_cast<Value>(0);
        for (const auto& rThreadStorage : loop.getPool().getStorage())
            output += std::get<0>(rThreadStorage.values());
        return output;
    } /*if _maybeThreads.has_value()*/ else {
        return std::inner_product(
            left.begin(),
            left.end(),
            right.begin(),
            static_cast<Value>(0));
    } // else
}


template <class T, TagLike TTag>
requires (std::is_same_v<TTag,tags::Serial> || std::is_same_v<TTag,tags::SMP>)
void DefaultSpace<T,TTag>::scale(VectorView target, ConstVectorView source, Value scale) const {
    const auto job = [target, source, this] (const auto& op) -> void {
        if (_maybeThreads.has_value()) {
            mp::ParallelFor<>(_maybeThreads.value())
                .execute(
                    this->size(target),
                    [target, source, &op] (std::size_t iComponent) -> void {
                        target[iComponent] = op(target[iComponent], source[iComponent]);
                    });
        } else {
            std::transform(
                target.begin(),
                target.end(),
                source.begin(),
                target.begin(),
                op);
        }
    }; // job

    if (scale == static_cast<T>(1))
        job(std::multiplies<T>());
    else if (scale == static_cast<T>(-1))
        job([] (Value left, Value right) -> Value {return -left * right;});
    else
        job([scale] (Value left, Value right) -> Value {return scale * left * right;});
}


template <class T, TagLike TTag>
requires (std::is_same_v<TTag,tags::Serial> || std::is_same_v<TTag,tags::SMP>)
void DefaultSpace<T,TTag>::scale(VectorView view, Value value) const {
    if (_maybeThreads.has_value()) {
        mp::ParallelFor<>(_maybeThreads.value())
            .execute(
                view,
                [value] (Ref<Value> rComponent) -> void {
                    rComponent *= value;
                });
    } else {
        for (Ref<Value> rComponent : view)
            rComponent *= value;
    }
}


template <class T, TagLike TTag>
requires (std::is_same_v<TTag,tags::Serial> || std::is_same_v<TTag,tags::SMP>)
void DefaultSpace<T,TTag>::add(VectorView target, ConstVectorView source, Value scale) const {
    const auto job = [&target, &source, this] (auto op) -> void {
        if (_maybeThreads.has_value()) {
            mp::ParallelFor<>(_maybeThreads.value())
                .execute(
                    target.size(),
                    [&op, &target, &source] (std::size_t iComponent) -> void {
                        target[iComponent] = op(target[iComponent], source[iComponent]);
                    });
        } else {
            std::transform(
                target.begin(),
                target.end(),
                source.begin(),
                target.begin(),
                op);
        }
    }; // job

    if (scale == static_cast<Value>(1))
        job([] (Value left, Value right) -> Value {return left + right;});
    else if (scale == static_cast<Value>(-1))
        job([] (Value left, Value right) -> Value {return left - right;});
    else
        job([scale] (Value left, Value right) -> Value {return left + scale * right;});
}


template <class T, TagLike TTag>
requires (std::is_same_v<TTag,tags::Serial> || std::is_same_v<TTag,tags::SMP>)
void DefaultSpace<T,TTag>::assign(VectorView target, ConstVectorView source) const {
    if (_maybeThreads.has_value()) {
        mp::ParallelFor<>(_maybeThreads.value())
            .execute(
                target.size(),
                [&target, &source] (std::size_t iComponent) -> void {
                    target[iComponent] = source[iComponent];
                });
    } else {
        std::copy(
            source.begin(),
            source.end(),
            target.begin());
    } // else
}


template <class T, TagLike TTag>
requires (std::is_same_v<TTag,tags::Serial> || std::is_same_v<TTag,tags::SMP>)
void DefaultSpace<T,TTag>::fill(VectorView view, Value value) const {
    if (_maybeThreads.has_value()) {
        mp::ParallelFor<>(_maybeThreads.value())
            .execute(
                view,
                [value] (Ref<Value> rComponent) -> void {
                    rComponent = value;
                });
    } else {
        std::fill(
            view.begin(),
            view.end(),
            value);
    } // else
}


template <class T, TagLike TTag>
requires (std::is_same_v<TTag,tags::Serial> || std::is_same_v<TTag,tags::SMP>)
OptionalRef<mp::ThreadPoolBase> DefaultSpace<T,TTag>::getThreads() {
    return _maybeThreads;
}


template class DefaultSpace<float,tags::Serial>;
template class DefaultSpace<float,tags::SMP>;
template class DefaultSpace<double,tags::Serial>;
template class DefaultSpace<double,tags::SMP>;


} // namespace cie::linalg
