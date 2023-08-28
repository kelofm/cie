// attributes: position,vertexColor
// sizes: 3,3
// strides: 6,6
// offsets: 0,3
// uniforms: transformation
// textures: 

#version 450 core

in vec3 position;
in vec3 vertexColor;
out vec3 vColor;

uniform mat4 transformation;

void main()
{
    vColor      = vertexColor;
    gl_Position = transformation * vec4( position, 1.0 );
}