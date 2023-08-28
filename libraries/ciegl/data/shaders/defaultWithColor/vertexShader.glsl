#version 450 core

// --- Attribute Declarations ---
in vec3 position;
in vec3 vertexColor;

// --- Output Declaratios ---
out vec3 vColor;

// --- Uniform Declarations ---
uniform mat4 transformation;


void main()
{
    vColor      = vertexColor;
    gl_Position = transformation * vec4( position, 1.0 );
}