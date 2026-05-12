// --- Internal Includes ---
#include "embeddedPoisson2D/configuration.hpp"
#include "embeddedPoisson2D/definitions.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"


namespace cie::fem {


void makeSchema(Ref<cie::io::JSONSchema> rSchema) {
    CIE_BEGIN_EXCEPTION_TRACING
        cie::io::JSONObject json(std::string {R"({
            "$schema" : "http://json-schema.org/draft-07/schema",
            "title" : "embedded-poisson-2d",
            "description" : "",
            "type" : "object",
            "properties" : {
                "dirichlet-1d" : {"$ref" : "/cie/fem/dirichlet-condition-1d"},
                "domain" : {"$ref" : "/cie/fem/embedded-domain"},
                "discretization" : {"$ref" : "/cie/fem/discretization"},
                "linear-system" : {"$ref" : "/cie/fem/linear-system"}
            },
            "required" : ["domain"],
            "additionalProperties" : false
        })"});
        rSchema = cie::io::JSONSchema(std::move(json));
    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::fem
