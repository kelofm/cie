// --- Internal Includes ---
#include "poisson2D/configuration.hpp"
#include "poisson2D/definitions.hpp"

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
                "domains" : {
                    "type" : "array",
                    "items" : {
                        "$ref" : "/cie/fem/embedded-domain"
                    }
                },
                "dirichlet-1d" : {"$ref" : "/cie/fem/dirichlet-condition-1d"},
                "discretization" : {"$ref" : "/cie/fem/discretization"},
                "linear-system" : {"$ref" : "/cie/fem/constrained-linear-system"},
                "write-schema" : {
                    "type" : ["string", "null"],
                    "default" : null,
                    "description" : "Path to write the configuration schema to, or null."
                },
                "write-applied-configuration" : {
                    "type" : ["string", "null"],
                    "default" : null,
                    "description" : "Path to write the applied configuration JSON to, or null."
                }
            },
            "default" : {
                "dirichlet-1d" : {},
                "domains" : [{"type" : "default"}],
                "discretization" : {},
                "linear-system" : {}
            },
            "required" : ["domains"],
            "additionalProperties" : false
        })"});
        rSchema = cie::io::JSONSchema(std::move(json));
    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::fem
