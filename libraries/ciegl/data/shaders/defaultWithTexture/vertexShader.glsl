#version 450 core

// --- Attribute Declarations ---
in vec3 position;
in vec2 textureCoordinates;

// --- Output Declaratios ---
out vec2 vTextureCoordinates;

// --- Uniform Declarations ---
uniform mat4 transformation;


void main()
{
    vTextureCoordinates = textureCoordinates;
    gl_Position = transformation * vec4( position, 1.0 );
}