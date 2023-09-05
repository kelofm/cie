#version 450 core

// --- Layout Definition ---
layout(triangles) in;
layout(triangle_strip, max_vertices=3 ) out;

// --- Input Declarations ---
in vec3 vPosition[];
in vec4 vColor[];

// --- Output Declarations ---
out vec3 fPosition;
out vec3 fNormal;
out vec4 fColor;


void main()
{
    vec3 normal = normalize(cross(
        vPosition[1] - vPosition[0],
        vPosition[2] - vPosition[0]
    ));

    for ( int i=0; i<3; ++i ) {
        fNormal     = normal;
        fColor      = vColor[i];
        fPosition   = vPosition[i];
        gl_Position = gl_in[i].gl_Position;
        EmitVertex();
    }

    EndPrimitive();
}