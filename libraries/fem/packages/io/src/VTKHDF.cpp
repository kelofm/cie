#ifdef CIE_ENABLE_HDF5

// --- External Includes ---
#include "H5Cpp.h"

// --- FEM Includes ---
#include "packages/io/inc/VTKHDF.hpp"

// --- STL Includes ---
#include <array>
#include <variant>


namespace cie::io {


template <class T>
concept H5EntityLike
=  std::is_same_v<T,H5::Group>
|| std::is_same_v<T,H5::DataSet>;


struct VTKHDF::Output::Impl {
    using Prefix = VTKHDF::Output::Prefix;


    Impl(Ref<const std::filesystem::path> rPath)
        : _file() {
            CIE_BEGIN_EXCEPTION_TRACING
            H5std_string h5Path;
            const auto& rPathString = rPath.string();
            h5Path.resize(rPathString.size());

            std::copy(
                rPathString.begin(),
                rPathString.end(),
                h5Path.begin());

            // Create or wipe the output file.
            _file = H5::H5File(
                h5Path,
                H5F_ACC_TRUNC);
            CIE_END_EXCEPTION_TRACING

            CIE_BEGIN_EXCEPTION_TRACING
            H5::Group rootGroup(_file.createGroup("/VTKHDF"));
            this->writeAttribute(
                rootGroup,
                "Type",
                "MultiBlockDataSet");
            const std::array<int,2> version {2, 4};
            this->writeAttribute(
                rootGroup,
                "Version",
                std::span<const int>(version));
            CIE_END_EXCEPTION_TRACING
    }


    H5::Group findGroup(Ref<const Prefix> rPrefix) {
        CIE_BEGIN_EXCEPTION_TRACING
            H5::Group group = _file.openGroup("/");
            std::string groupName;

            for (const auto& rGroupName : rPrefix) {
                const auto& rString = rGroupName.string();
                groupName.resize(rString.size());
                std::copy(
                    rString.begin(),
                    rString.end(),
                    groupName.begin());
                group = group.openGroup(groupName.c_str());
            } // for rGroupName in rPrefix

            return group;
        CIE_END_EXCEPTION_TRACING
    }


    H5::DataSet findDataset(Ref<const Prefix> rPrefix) {
        CIE_BEGIN_EXCEPTION_TRACING
            H5::Group parent = this->findGroup(rPrefix.parent_path());
            const auto& rDatasetName = rPrefix.filename().string();
            std::string datasetString;
            datasetString.resize(rDatasetName.size());
            std::copy(
                rDatasetName.begin(),
                rDatasetName.end(),
                datasetString.begin());
            return parent.openDataSet(datasetString.c_str());
        CIE_END_EXCEPTION_TRACING
    }


    std::variant<H5::Group,H5::DataSet> find(Ref<const Prefix> rPrefix) {
        CIE_BEGIN_EXCEPTION_TRACING
            H5::Group parent = this->findGroup(rPrefix.parent_path());
            const auto& rChildName = rPrefix.filename().string();
            std::string childName;
            childName.resize(rChildName.size());
            std::copy(
                rChildName.begin(),
                rChildName.end(),
                childName.begin());

            const auto childType = parent.childObjType(childName.c_str());
            if (childType == H5O_TYPE_GROUP) {
                return parent.openGroup(childName.c_str());
            } else if (childType == H5O_TYPE_DATASET) {
                return parent.openDataSet(childName.c_str());
            } else CIE_THROW(Exception, rPrefix << " is of invalid type")
        CIE_END_EXCEPTION_TRACING
    }


    template <class T>
    static H5::PredType getH5Type() {
        if constexpr (std::is_same_v<T,int>) {
            if constexpr (sizeof(int) == 4) return H5::PredType::NATIVE_INT32;
            else if constexpr (sizeof(int) == 8) return H5::PredType::NATIVE_INT64;
            else static_assert(std::is_same_v<T,void>, "unsupported type");
        } else if constexpr (std::is_same_v<T,unsigned>) {
            if constexpr (sizeof(unsigned) == 4) return H5::PredType::NATIVE_UINT32;
            else if constexpr (sizeof(unsigned) == 8) return H5::PredType::NATIVE_UINT64;
            else static_assert(std::is_same_v<T,void>, "unsupported type");
        } else if constexpr (std::is_same_v<T,std::size_t>) {
            if constexpr (sizeof(std::size_t) == 4) return H5::PredType::NATIVE_UINT32;
            else if constexpr (sizeof(std::size_t) == 8) return H5::PredType::NATIVE_UINT64;
            else static_assert(std::is_same_v<T,void>, "unsupported type");
        } else if constexpr (std::is_same_v<T,std::int8_t>) return H5::PredType::NATIVE_INT8;
        else if constexpr (std::is_same_v<T,float>) return H5::PredType::IEEE_F32LE;
        else if constexpr (std::is_same_v<T,double>) return H5::PredType::IEEE_F64LE;
        else static_assert(std::is_same_v<T,void>, "unsupported_type");
    }


    template <H5EntityLike TEntity>
    static void writeAttribute(
        Ref<TEntity> rEntity,
        Ref<const std::string> rName,
        std::string_view value) {
            CIE_BEGIN_EXCEPTION_TRACING
                H5::StrType stringType(0, value.size());
                H5::DataSpace dataSpace(H5S_SCALAR);
                H5::Attribute attribute = rEntity.createAttribute(
                    rName.c_str(),
                    stringType,
                    dataSpace);
                attribute.write(
                    stringType,
                    value.data());
            CIE_END_EXCEPTION_TRACING
    }


    template <H5EntityLike TEntity, class TValue>
    static void writeAttribute(
        Ref<TEntity> rEntity,
        Ref<const std::string> rName,
        std::span<const TValue> value) {
            CIE_BEGIN_EXCEPTION_TRACING
                const auto h5Type = Impl::getH5Type<TValue>();
                std::array<hsize_t,1> shape {value.size()};
                H5::DataSpace dataSpace(
                    shape.size(),
                    shape.data());
                H5::Attribute attribute = rEntity.createAttribute(
                    rName.c_str(),
                    h5Type,
                    dataSpace);
                attribute.write(
                    h5Type,
                    value.data());
            CIE_END_EXCEPTION_TRACING
    }

private:
    H5::H5File _file;
}; // struct Impl


VTKHDF::Output::Output()
    : Output("output.vtkhdf")
{}


VTKHDF::Output::Output(Ref<const std::filesystem::path> rPath)
    : _pImpl(new Impl(rPath))
{}


VTKHDF::Output::~Output() = default;


void VTKHDF::Output::writeAttribute(
    Ref<const Prefix> rPrefix,
    Ref<const std::string> rName,
    std::string_view value) {
        CIE_BEGIN_EXCEPTION_TRACING
        std::visit(
            [&] (auto entity) -> void {
                _pImpl->writeAttribute(
                    entity,
                    rName,
                    value);
            },
            _pImpl->find(rPrefix));
        CIE_END_EXCEPTION_TRACING
}


template <class T>
void VTKHDF::Output::writeAttribute(
    Ref<const Prefix> rPrefix,
    Ref<const std::string> rName,
    std::span<const T> value) {
        CIE_BEGIN_EXCEPTION_TRACING
        std::visit(
            [&] (auto entity) -> void {
                _pImpl->writeAttribute(
                    entity,
                    rName,
                    value);
            },
            _pImpl->find(rPrefix));
        CIE_END_EXCEPTION_TRACING
}


#define CIE_INSTANTIATE_VTKHDF(T)                       \
    template void VTKHDF::Output::writeAttribute<T>(    \
        Ref<const Prefix>,                              \
        Ref<const std::string>,                         \
        std::span<const T>);


CIE_INSTANTIATE_VTKHDF(int)
CIE_INSTANTIATE_VTKHDF(unsigned)
CIE_INSTANTIATE_VTKHDF(std::size_t)
CIE_INSTANTIATE_VTKHDF(float)
CIE_INSTANTIATE_VTKHDF(double)


#undef CIE_INSTANTIATE_VTKHDF


} // namespace cie::io


#endif
