#version 450 core

// --- InputDeclarations ---
in vec3 fColor;

// --- Output Declarations ---
out vec4 color;


void main()
{
    color = vec4( fColor, 1.0 );
}