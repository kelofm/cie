defaultVertexShaderWithTexture = '''

varying vec3 position;
varying vec3 normal;
varying vec2 textureCoordinates;


void main() {
    normal              = $normal;
    textureCoordinates  = $textureCoordinates;
    position            = $position;
    gl_Position         = $transform(vec4(position,1));
}

'''


incompleteAggregatedVertexShader = '''
#version 130

// Generated declarations ----------------------------------------------------------------
uint textureIDMap[<textureSize_placeholder>] = uint[<textureSize_placeholder>](<textureIndex_placeholder>);

uint textureOffsetX[<textureSize_placeholder>] = uint[<textureSize_placeholder>](<textureOffsetX_placeholder>);
uint textureOffsetY[<textureSize_placeholder>] = uint[<textureSize_placeholder>](<textureOffsetY_placeholder>);
uint textureWidth[<textureSize_placeholder>] = uint[<textureSize_placeholder>](<textureWidth_placeholder>);
uint textureHeight[<textureSize_placeholder>] = uint[<textureSize_placeholder>](<textureHeight_placeholder>);


// Static declarations -------------------------------------------------------------------
varying vec3 position;
varying vec3 normal;
varying vec2 textureCoordinates;

uint aggregatedTextureWidth     = uint($aggregatedTextureWidth);
uint aggregatedTextureHeight    = uint($aggregatedTextureHeight);

varying float materialID;


// Function definitions ------------------------------------------------------------------
vec2 getGlobalTextureCoordinates( uint matID, vec2 localTextureCoordinates )
{
    uint index  = textureIDMap[matID];
    return vec2(    ( float(textureOffsetY[index]) + float(textureHeight[index]) * (1.0-localTextureCoordinates.y) ) / float(aggregatedTextureHeight),
                    ( float(textureOffsetX[index]) + float(textureWidth[index]) * localTextureCoordinates.x ) / float(aggregatedTextureWidth) );
}


// MAIN ----------------------------------------------------------------------------------
void main() {
    materialID          = $materialID;
    normal              = $normal;
    textureCoordinates  = getGlobalTextureCoordinates( uint(materialID), vec2($textureCoordinates) );
    position            = $position;
    gl_Position         = $transform(vec4(position,1));
}

'''