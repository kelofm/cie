#pragma once

// --- Utility Includes ---
#include "packages/concurrency/inc/ThreadPoolBase.hpp"
#include "packages/stl_extension/inc/OptionalRef.hpp"

// --- Linalg Includes ---
#include "packages/utilities/inc/CSRView.hpp"


namespace cie::linalg {


template <class TValue, class TIndex>
void CSRSum(
    CSRView<TValue,TIndex> left,
    CSRView<const TValue,const TIndex> right,
    TValue scale = static_cast<TValue>(1),
    OptionalRef<mp::ThreadPoolBase> rMaybeThreads = {});


} // namespace cie::linalg
