#ifndef CIE_GL_GLM_OVERLOADS_IMPL_HPP
#define CIE_GL_GLM_OVERLOADS_IMPL_HPP

// --- STL Includes ---
#include <sstream>


namespace cie::utils {


template <typename T, glm::qualifier Q>
Loggee& operator<<( Loggee& r_loggee,
                    const glm::vec<3,T,Q>& r_vector )
{
    std::stringstream stream;
    for ( Size i=0; i<Size(r_vector.length()); ++i )
        stream << r_vector[i] << ' ';
    
    return r_loggee << stream.str();
}


} // namespace cie::utils

#endif