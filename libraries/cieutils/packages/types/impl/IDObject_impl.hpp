#ifndef CIE_UTILS_ID_OBJECT_IMPL_HPP
#define CIE_UTILS_ID_OBJECT_IMPL_HPP


namespace cie::utils {


template <class TID>
IDObject<TID>::IDObject() :
    IDObject<TID>(0)
{
}


template <class TID>
IDObject<TID>::IDObject(const TID& r_id) :
    _id(r_id)
{
}


template <class TID>
inline void
IDObject<TID>::setID(TID id)
{
    _id = id;
}


template <class TID>
inline const TID
IDObject<TID>::getID() const
{
    return _id;
}


} // namespace cie::utils


#endif