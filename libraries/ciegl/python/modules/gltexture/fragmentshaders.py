defaultFragmentShaderWithTexture = '''

varying vec3 position;
varying vec3 normal;
varying vec2 textureCoordinates;

uniform sampler2D u_objectTexture;


void main() {
    // Texture
    vec4 color                  = texture2D( u_objectTexture, textureCoordinates );

    // Diffuse
    vec3 lightDir               = normalize( vec3($lightPos) - position );
    float diffuse               = dot( normal, lightDir );
    diffuse                     = min( max(diffuse,0.0), 1.0 );

    // Specular
    vec3 halfWayVector          = normalize(  ( lightDir + normalize(position - vec3($cameraPos)) )/2.0  );
    float specular              = min( max(   dot(halfWayVector,normal)   ,0.0), 1.0 );
    specular                    = specular * specular * specular * specular;
    specular                    = specular * specular * specular * specular;
    specular                    = specular * specular * specular * specular;

    // Additive
    vec3 fragColor              = vec3($lightColor) * color.xyz;
    fragColor                   = fragColor * ( $diffuseMaterialConstant * diffuse +
                                                $specularMaterialConstant * specular)
                                                + $ambientMaterialConstant * vec3($ambientLight) * color.xyz;


    gl_FragColor                = vec4(fragColor,1.0);
}

'''


incompleteAggregatedFragmentShader    = '''
#version 130

// Generated declarations ----------------------------------------------------------------
uint textureIDMap[<textureSize_placeholder>] = uint[<textureSize_placeholder>](<textureIndex_placeholder>);

float ambientConstants[<textureSize_placeholder>] = float[<textureSize_placeholder>](<ambientConstants_placeholder>);
float diffuseConstants[<textureSize_placeholder>] = float[<textureSize_placeholder>](<diffuseConstants_placeholder>);
float specularConstants[<textureSize_placeholder>] = float[<textureSize_placeholder>](<specularConstants_placeholder>);


// Static declarations -------------------------------------------------------------------
varying vec3 position;
varying vec3 normal;
varying vec2 textureCoordinates;

uniform sampler2D aggregatedTexture;

varying float materialID;


// MAIN ----------------------------------------------------------------------------------
void main() {

    // Material index
    uint materialIndex          = textureIDMap[int(materialID)];

    // Texture
    vec4 color                  = texture2D( aggregatedTexture, textureCoordinates );

    // Ambient
    vec3 ambient                = ambientConstants[materialIndex] * vec3($ambientLight) * color.xyz;

    // Diffuse
    vec3 lightDir               = normalize( vec3($lightPos) - position );
    float diffuse               = dot( normal, lightDir );
    diffuse                     = min( max(diffuse,0.0), 1.0 );
    diffuse                     = diffuse * diffuseConstants[materialIndex];

    // Specular
    vec3 halfWayVector          = normalize(  ( lightDir + normalize(position - vec3($cameraPos)) )/2.0  );
    float specular              = min( max(   dot(halfWayVector,normal)   ,0.0), 1.0 );
    specular                    = specular * specular * specular * specular;
    specular                    = specular * specular * specular * specular;
    specular                    = specular * specular * specular * specular;
    specular                    = specular * specularConstants[materialIndex];

    // Additive
    vec3 fragColor              = vec3($lightColor) * color.xyz;
    fragColor                   = fragColor * (diffuse + specular) + ambient;

    gl_FragColor                = vec4(fragColor,1.0);
}


'''