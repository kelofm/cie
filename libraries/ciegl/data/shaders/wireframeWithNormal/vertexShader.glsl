#version 450 core

// --- Input Declarations ---
in vec3 position;
in vec3 normal;

// --- Uniform Declarations ---
uniform mat4 transformation;
uniform float wireframeScale;


void main()
{
    gl_Position = transformation * vec4(position + wireframeScale * normal, 1.0);
}