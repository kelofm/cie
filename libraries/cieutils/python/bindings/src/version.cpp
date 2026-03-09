// --- External Includes ---
#include "pybind11/pybind11.h"
#include "pybind11/operators.h"

// --- Utility Includes ---
#include "version.hpp"
#include "python/bindings/inc/common.hpp"

// --- STL Includes ---
#include <sstream>


namespace cie {


template <>
void addBindings<Version>(Ref<pybind11::module_> rModule) {
    pybind11::class_<cie::Version>(rModule, "Version")
        .def(pybind11::init<cie::Ref <const cie::Version::Branch>, cie::Ref<const cie::Version::Hash>>())
        .def("__str__", [](cie::Ref<const cie::Version> rVersion){return (std::string) rVersion;})
        .def(pybind11::self == pybind11::self)
        .def(pybind11::self != pybind11::self)
        .def(pybind11::self < pybind11::self)
        .def_property_static("local", [](pybind11::object) {return cie::Version::local;}, [](cie::Version) {})
        ;
}


} // namespace cie
