#version 450 core

// --- Layout Definition ---
layout(triangles) in;
layout(line_strip, max_vertices=4 ) out;

// --- Input Definition ---
in vec3 vNormal[];

// --- Output Definition ---
out vec3 fNormal;


void main()
{
    for( int i=0; i<4; ++i )
    {
        fNormal     = vNormal[i%3];
        gl_Position = gl_in[i%3].gl_Position;
        EmitVertex();
    }

    EndPrimitive();
}