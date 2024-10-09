#ifndef CIE_UTILS_ID_OBJECT_HPP
#define CIE_UTILS_ID_OBJECT_HPP

// --- Internal Includes ---
#include "packages/types/inc/types.hpp"


namespace cie::utils {


/// @ingroup cieutils
template <class TID = Size>
class IDObject
{
public:
    using id_type = TID;

public:
    IDObject();

    IDObject(const TID& r_id);

    IDObject(IDObject&& r_rhs) = default;

    IDObject(const IDObject& r_rhs) = default;

    IDObject& operator=(IDObject&& r_rhs) = default;

    IDObject& operator=(const IDObject& rhs) = default;

    void setID(TID id);

    const TID getID() const;

protected:
    TID _id;
};


} // namespace cie::utils

#include "packages/types/impl/IDObject_impl.hpp"

#endif