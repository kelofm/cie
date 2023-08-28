#version 450 core

// --- Layout Definition ---
layout(triangles) in;
layout(triangle_strip, max_vertices=3 ) out;

// --- Input Declarations ---
in vec3 vPosition[];
in vec3 vNormal[];

// --- Output Declarations ---
out vec3 fPosition;
out vec3 fNormal;


void main()
{
    for ( int i=0; i<3; ++i )
    {
        fNormal   = vNormal[i];
        fPosition = vPosition[i];
        gl_Position = gl_in[i].gl_Position;
        EmitVertex();
    }

    EndPrimitive();
}