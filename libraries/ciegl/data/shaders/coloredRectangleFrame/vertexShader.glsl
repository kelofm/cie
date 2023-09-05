#version 450 core

// --- Input Declarations ---
in vec2 position;
in vec2 lengths;
in vec4 color;

// --- Output Declarations ---
out vec2 vLengths;
out vec4 vColor;


void main()
{
    vLengths = lengths;
    vColor = color;
    gl_Position = vec4(position, 0.0, 1.0);
}
