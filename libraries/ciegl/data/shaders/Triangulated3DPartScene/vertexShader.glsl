#version 450 core

// --- Input Declarations ---
in vec3 position;

// --- Uniform Declarations ---
uniform mat4 transformation;

// --- Output Declarations ---
out vec3 vPosition;


void main()
{
    vPosition = position;
    gl_Position = vec4( transformation * vec4(position, 1.0) );
}