
const drawTimeElem = document.querySelector("#drawTime");
var canvas = document.getElementById('canvas');

gl = canvas.getContext('webgl2');
if (!gl)
{
	console.error("webgl2 not available");
}

var vertices =
[
	-1.0, 2.0, 0.0,
	-1.0, -1.0, 0.0,
	2.0, -1.0, 0.0,
];

indices = [0, 1, 2];

var vertexBuffer = gl.createBuffer();
gl.bindBuffer(gl.ARRAY_BUFFER, vertexBuffer);
gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(vertices), gl.STATIC_DRAW);
gl.bindBuffer(gl.ARRAY_BUFFER, null);

var indexBuffer = gl.createBuffer();
gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, indexBuffer);
gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(indices), gl.STATIC_DRAW);
gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, null);

var vertexShaderCode = `
precision highp float;

attribute vec3 coordinates;

void main(void)
{
	gl_Position = vec4(coordinates, 0.000001);
}
`;

var vertexShader = gl.createShader(gl.VERTEX_SHADER);
gl.shaderSource(vertexShader, vertexShaderCode);
gl.compileShader(vertexShader);

var fragmentShaderCode = `
precision highp float;

uniform float width;
uniform float height;
uniform float time;

void main(void)
{
	gl_FragColor = vec4(sin(time / 1000.0) * gl_FragCoord.x / width, gl_FragCoord.y / height, gl_FragCoord.z, 1.0);
}
`;

var fragmentShader = gl.createShader(gl.FRAGMENT_SHADER);
gl.shaderSource(fragmentShader, fragmentShaderCode);
gl.compileShader(fragmentShader);

var shaderProgram = gl.createProgram();
gl.attachShader(shaderProgram, vertexShader);
gl.attachShader(shaderProgram, fragmentShader);
gl.linkProgram(shaderProgram);
gl.useProgram(shaderProgram);

gl.bindBuffer(gl.ARRAY_BUFFER, vertexBuffer);
gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, indexBuffer);

var coord = gl.getAttribLocation(shaderProgram, "coordinates");
gl.vertexAttribPointer(coord, 3, gl.FLOAT, false, 0, 0);
gl.enableVertexAttribArray(coord);

gl.clearColor(0.0, 0.0, 0.0, 1.0);
gl.enable(gl.DEPTH_TEST);
gl.clear(gl.COLOR_BUFFER_BIT);
gl.viewport(0, 0, canvas.width, canvas.height);

var widthHandle = gl.getUniformLocation(shaderProgram, "width");
var heightHandle = gl.getUniformLocation(shaderProgram, "height");
var timeHandle = gl.getUniformLocation(shaderProgram, "time");

gl.uniform1f(widthHandle, canvas.width);
gl.uniform1f(heightHandle, canvas.height);
gl.uniform1f(timeHandle, 0.0);

function changeShader(now)
{
	gl.detachShader(shaderProgram, fragmentShader);
	var fragmentShaderCode = `
precision highp float;

uniform float width;
uniform float height;
uniform float time;

void main(void)
{
	gl_FragColor = vec4(abs(sin(time / 1000.0)) * gl_FragCoord.x / width, gl_FragCoord.y / height, gl_FragCoord.z, 1.0);
}
`;
	gl.shaderSource(fragmentShader, fragmentShaderCode);
	gl.compileShader(fragmentShader);
	gl.attachShader(shaderProgram, fragmentShader);
	gl.linkProgram(shaderProgram);
	gl.useProgram(shaderProgram);
	widthHandle = gl.getUniformLocation(shaderProgram, "width");
	heightHandle = gl.getUniformLocation(shaderProgram, "height");
	timeHandle = gl.getUniformLocation(shaderProgram, "time");
	gl.uniform1f(widthHandle, canvas.width);
	gl.uniform1f(heightHandle, canvas.height);
	gl.uniform1f(timeHandle, now);
}

var then = 0;
var delta = 0;
function drawScene(now)
{
	gl.drawElements(gl.TRIANGLES, indices.length, gl.UNSIGNED_SHORT, 0);
	changeShader(now);
	delta = now - then;
	then = now;
	drawTimeElem.textContent = delta.toFixed(3);
	requestAnimationFrame(drawScene);
}

requestAnimationFrame(drawScene);
