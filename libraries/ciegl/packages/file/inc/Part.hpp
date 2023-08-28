#ifndef CIE_GL_PART_HPP
#define CIE_GL_PART_HPP

// --- Utility Includes ---
#include "packages/types/inc/types.hpp"

// --- Internal Includes ---
#include "packages/buffer/inc/Buffer.hpp"
#include "packages/traits/inc/GLTraits.hpp"

// --- STL Includes ---
#include <vector>
#include <memory>
#include <istream>
#include <filesystem>


namespace cie::gl {


class Part : public GLTraits
{
public:
    using GLTraits::attribute_type;
    using GLTraits::attribute_container;

    using GLTraits::index_type;
    using GLTraits::index_container;

public:
    Part();

    /// Total size of all attributes in bytes
    Size attributeByteCount() const;

    /// Total size of all indices in bytes
    Size indexByteCount() const;

    /// Get number of stored vertices
    Size numberOfVertices() const;

    /// Get number of stored primitives
    Size numberOfPrimitives() const;

    /// Load data from input stream
    virtual void load( std::istream& r_stream );

    /// Get attribute container
    const attribute_container::SharedPointer& attributes() const;

    /// Get primitive vertex maps
    const index_container& indices() const;

    /// Get number of dimensions
    virtual Size dimension() const = 0;

    /// Get number of attributes that make up a vertex
    virtual Size vertexAttributeSize() const = 0;

    /// Get number of vertices that make up a primitive
    virtual Size primitiveVertexSize() const = 0;

    /// Get number of attributes that make up a primitive
    Size primitiveAttributeSize() const;

    /// Get number of bytes needed for defining a single primitive
    Size primitiveByteSize() const;

    /// Write part in STL format
    virtual void writeSTL( const std::filesystem::path& r_path ) const = 0;

protected:
    Part( const Part& r_rhs ) = delete;
    Part( Part&& r_rhs ) = delete;
    Part& operator=( const Part& r_rhs ) = delete;

protected:
    attribute_container::SharedPointer _p_attributes;
    index_container _indices;
};


using PartPtr = std::shared_ptr<Part>;


} // namespace cie::gl


#endif