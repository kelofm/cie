#version 450 core

// --- Input Declarations ---
in vec3 position;

// --- Uniform Declarations ---
uniform mat4 transformation;


void main()
{
    gl_Position = transformation * vec4( position, 1.0 );
}
