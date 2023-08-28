// attributes: position,cellSize
// sizes: 3,1
// strides: 4,4
// offsets: 0,3
// uniforms: transformation
// textures: 

#version 450 core

// Vertex layout:
//  poition      : x, y, z
//  cell size    : size


in vec3     position;
in float    cellSize;
out float   vCellSize;

uniform mat4 transformation;

void main()
{
    vCellSize   = cellSize;
    gl_Position = transformation * vec4( position, 1.0 );
}