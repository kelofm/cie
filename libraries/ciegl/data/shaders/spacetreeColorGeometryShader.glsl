// attributes: 
// sizes: 
// strides: 
// offsets: 
// uniforms: 
// textures: 

#version 450 core

layout(triangles) in;
layout(triangle_strip, max_vertices=3 ) out;

in vec3 vColor[];
out vec3 fColor;

void main()
{
    fColor      = vColor[0];
    gl_Position = gl_in[0].gl_Position;
    EmitVertex();

    fColor      = vColor[1];
    gl_Position = gl_in[1].gl_Position;
    EmitVertex();

    fColor      = vColor[2];
    gl_Position = gl_in[2].gl_Position;
    EmitVertex();

    EndPrimitive();
}