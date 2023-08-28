#version 450 core

// --- Layout Definition ---
layout(points) in;
layout(line_strip, max_vertices=5 ) out;

// --- Uniform Definitions ---
uniform mat4 transformation;

// --- Input Definition ---
in vec2[] vLengths;


void main()
{
    gl_Position = transformation * gl_in[0].gl_Position;
    EmitVertex();

    gl_Position = transformation * (gl_in[0].gl_Position + vec4( vLengths[0][0], 0.0, 0.0, 0.0 ) );
    EmitVertex();

    gl_Position = transformation * (gl_in[0].gl_Position + vec4( vLengths[0], 0.0, 0.0 ) );
    EmitVertex();

    gl_Position = transformation * (gl_in[0].gl_Position + vec4( 0.0, vLengths[0][1], 0.0, 0.0 ) );
    EmitVertex();

    gl_Position = transformation * gl_in[0].gl_Position;
    EmitVertex();

    EndPrimitive();
}