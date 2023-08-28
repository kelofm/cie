defaultVertexShader = '''

varying vec3 position;
varying vec3 normal;
varying vec4 color;


void main() {
    normal      = $normal;
    color       = $color;
    position    = $position;
    gl_Position = $transform(vec4(position,1));
}

'''