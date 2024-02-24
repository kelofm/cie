// --- FEM Includes ---
#include "packages/graph/inc/Assembler.hpp"


namespace cie::fem {


Assembler::Assembler() noexcept
    : _dofCounter(0)
{
}


Assembler::Assembler(std::size_t dofBegin) noexcept
    : _dofCounter(dofBegin)
{
}


std::size_t Assembler::dofCount() const noexcept
{
    return _dofCounter;
}


} // namespace cie::fem
