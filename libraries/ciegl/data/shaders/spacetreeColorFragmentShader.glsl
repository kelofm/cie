// attributes: color
// sizes: 4
// strides: 0
// offsets: 0
// uniforms: 
// textures: 

#version 450 core

in vec3 fColor;
out vec4 color;

void main()
{
    color = vec4( fColor, 1.0 );
}