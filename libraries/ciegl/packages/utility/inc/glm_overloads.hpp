#ifndef CIE_GL_GLM_OVERLOADS_HPP
#define CIE_GL_GLM_OVERLOADS_HPP

// --- External Includes ---
#include "glm/detail/type_vec3.hpp"

// --- Utility Includes ---
#include "packages/logging/inc/Loggee.hpp"


namespace cie::utils {


template <typename T, glm::qualifier Q>
Loggee& operator<<( Loggee& r_loggee,
                    const glm::vec<3,T,Q>& r_vector );


} // namespace cie::utils

#include "packages/utility/impl/glm_overloads.hpp"

#endif