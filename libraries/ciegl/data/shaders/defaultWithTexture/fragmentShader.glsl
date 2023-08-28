#version 450 core

// --- Input Declarations ---
in vec2 fTextureCoordinates;

// --- Uniform Declarations ---
uniform sampler2D image;

// --- Output Declarations ---
out vec4 color;


void main()
{
    color = texture(image, fTextureCoordinates);
}