#version 450 core

// --- Input Declarations ---
in vec3 position;
in vec3 normal;

// --- Uniform Declarations ---
uniform mat4 transformation;

// --- Output Declarations ---
out vec3 vPosition;
out vec3 vNormal;


void main()
{
    vPosition = position;
    vNormal   = normal;
    gl_Position = vec4( transformation * vec4(position, 1.0) );
}