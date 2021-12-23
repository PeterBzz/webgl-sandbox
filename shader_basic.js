
const shader_basic = new ShaderHelper(
// FRAGMENT SHADER
`#version 300 es

#ifdef GL_ES
precision highp float;
#endif

in vec4 color;

out vec4 fragColor;

void main(void)
{
	fragColor = vec4(0.5, 0.5, 0.5, 1.0);
}
`,
// VERTEX SHADER
`#version 300 es

in vec3 coordinates;

out vec4 color;

const vec4 white = vec4(0.5);

void main(void)
{
	color = white;
	gl_Position = vec4(coordinates, 0.000001);
}
`,
// CANVAS ID
'canvas',
// DRAW TIME ID
'drawTime'
);
