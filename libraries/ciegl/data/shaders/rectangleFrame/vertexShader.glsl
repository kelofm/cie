#version 450 core

// --- Input Declarations ---
in vec2 position;
in vec2 lengths;

// --- Output Declarations ---
out vec2 vLengths;


void main()
{
    vLengths = lengths;
    gl_Position = vec4(position, 0.0, 1.0);
}
