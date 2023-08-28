#version 450 core

// --- Input Declarations ---
in vec3 position;
in vec4 color;

// --- Uniform Declarations ---
uniform mat4 transformation;

// --- Output Declarations ---
out vec3 vPosition;
out vec4 vColor;


void main()
{
    vPosition = position;
    vColor    = color;
    gl_Position = vec4( transformation * vec4(position, 1.0) );
}