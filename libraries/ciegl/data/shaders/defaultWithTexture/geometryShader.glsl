#version 450 core

// --- Layout Definition ---
layout(triangles) in;
layout(triangle_strip, max_vertices=3 ) out;

// --- Input Declarations ---
in vec2 vTextureCoordinates[];

// --- Output Declarations ---
out vec2 fTextureCoordinates;


void main()
{
    for ( int i=0; i<3; ++i )
    {
        fTextureCoordinates = vTextureCoordinates[i];
        gl_Position = gl_in[i].gl_Position;
        EmitVertex();
    }

    EndPrimitive();
}