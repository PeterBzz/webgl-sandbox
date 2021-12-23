
const shader_colors = new ShaderHelper(
// FRAGMENT SHADER
`#version 300 es

#ifdef GL_ES
precision highp float;
#endif

uniform float width;
uniform float height;
uniform float time;

out vec4 fragColor;

void main(void)
{
	fragColor = vec4(sin(time / 1000.0) * gl_FragCoord.x / width, gl_FragCoord.y / height, gl_FragCoord.z, 1.0);
}
`,
// VERTEX SHADER
`#version 300 es

in vec3 coordinates;

void main(void)
{
	gl_Position = vec4(coordinates, 0.000001);
}
`,
// CANVAS ID
'canvas',
// DRAW TIME ID
'drawTime'
);
