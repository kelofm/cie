#version 450 core

// --- Input Declarations ---
in vec3 fPosition;
in vec3 fNormal;

// --- Uniform Declarations ---
uniform vec3 cameraPosition;

// --- Output Declarations
out vec4 color;


void main()
{
    float colorIntensity = dot(
        normalize( fNormal ),
        normalize( cameraPosition - fPosition )
    );

    float ambient = 0.1;

    colorIntensity = (1.0 - ambient) * max( 0.0, colorIntensity );
    colorIntensity *= colorIntensity * colorIntensity;
    colorIntensity *= colorIntensity * colorIntensity;
    colorIntensity += ambient;

    color = vec4( colorIntensity, colorIntensity, colorIntensity, 1.0 );
}