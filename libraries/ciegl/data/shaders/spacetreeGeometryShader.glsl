// attributes: 
// sizes: 
// strides: 
// offsets: 
// uniforms: 
// textures: 

#version 450 core

layout(points) in;
layout(line_strip, max_vertices=6 ) out;

in float vCellSize[];

void main()
{
    gl_Position = gl_in[0].gl_Position - vCellSize[0]/2.0f * vec4(1.0f, 0.0f, 0.0f, 1.0f);
    EmitVertex();
    gl_Position = gl_in[0].gl_Position + vCellSize[0]/2.0f * vec4(1.0f, 0.0f, 0.0f, 1.0f);
    EmitVertex();
    EndPrimitive();

    gl_Position = gl_in[0].gl_Position - vCellSize[0]/2.0f * vec4(0.0f, 1.0f, 0.0f, 1.0f);
    EmitVertex();
    gl_Position = gl_in[0].gl_Position + vCellSize[0]/2.0f * vec4(0.0f, 1.0f, 0.0f, 1.0f);
    EmitVertex();
    EndPrimitive();

    gl_Position = gl_in[0].gl_Position - vCellSize[0]/2.0f * vec4(0.0f, 0.0f, 1.0f, 1.0f);
    EmitVertex();
    gl_Position = gl_in[0].gl_Position + vCellSize[0]/2.0f * vec4(0.0f, 0.0f, 1.0f, 1.0f);
    EmitVertex();
    EndPrimitive();
}