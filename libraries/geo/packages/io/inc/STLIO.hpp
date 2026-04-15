#pragma once

// --- Utility Includes ---
#include "packages/types/inc/types.hpp"

// --- STL Includes ---
#include <span>
#include <memory>
#include <iosfwd>


namespace cie::io {


struct STLIO {
    template <class TValue, unsigned Dimension>
    class Input {
    public:
        Input() noexcept = default;

        Input(Ref<std::istream> rStream);

        ~Input();

        std::size_t triangleCount() const noexcept;

        void execute(
            std::span<TValue> coordinates,
            std::span<TValue> normals = {});

    private:
        struct Impl;
        std::unique_ptr<Impl> _pImpl;
    }; // class Input

    template <class TValue, unsigned Dimension>
    class Output {
    public:
        Output() noexcept = default;

        Output(Ref<std::ostream> rStream);

        ~Output();

        void execute(
            std::span<const TValue> coordinates,
            std::span<const TValue> normals = {});

    private:
        struct Impl;
        std::unique_ptr<Impl> _pImpl;
    }; // class Output
}; // struct STLIO


} // namespace cie::io
