#version 450 core

// --- Input Declarations ---
in vec2 position;

// --- Uniform Declarations ---
uniform mat4 transformation;


void main()
{
    gl_Position = transformation * vec4( position, 0.0, 1.0 );
}
