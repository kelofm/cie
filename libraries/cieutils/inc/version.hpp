#ifndef CIE_UTILS_VERSION_HPP

// --- Utility Includes ---
#include "packages/types/inc/types.hpp"

// --- STL Includes ---
#include <string>
#include <ostream>


namespace cie {


class Version
{
public:
    using Branch = std::string;

    using Hash = std::string;

public:
    static const Version local;

public:
    Version() = delete;

    Version(RightRef<Branch> r_branch, RightRef<Hash> r_hash) noexcept;

    Version(Ref<const Branch> r_branch, Ref<const Hash> r_hash);

    Version(Version&& r_rhs) noexcept = default;

    Version(const Version& r_rhs) = default;

    explicit operator std::string () const;

    friend Bool operator<(Ref<const Version> r_lhs, Ref<const Version> r_rhs);

    friend Bool operator==(Ref<const Version> r_lhs, Ref<const Version> r_rhs);

    friend Ref<std::ostream> operator<<(Ref<std::ostream> r_stream, Ref<const Version> r_version);

private:
    Branch _branch;

    Hash _hash;
}; // class Version


} // namespace cie


#endif
