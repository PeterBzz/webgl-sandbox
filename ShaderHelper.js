
class ShaderHelper
{
	static vertices =
	[
		-1.0, 2.0, 0.0,
		-1.0, -1.0, 0.0,
		2.0, -1.0, 0.0,
	];
	static indices = [0, 1, 2];
	static vertexBuffer = null;
	static indexBuffer = null;

	#drawTimeElement = null;
	#canvasElement = null;
	#context = null;
	#vertexShader = null;
	#fragmentShader = null;
	#program = null;
	#requestId = null;

	constructor(fragmentShaderCode, vertexShaderCode, canvasId, drawTimeId)
	{
		this.drawTimeElement = ShaderHelper.GetElement(drawTimeId);
		this.canvasElement = ShaderHelper.GetElement(canvasId);
		this.context = ShaderHelper.GetContext(this.canvasElement, 'webgl2');
		this.fragmentShaderCode = fragmentShaderCode;
		this.vertexShaderCode = vertexShaderCode;
		this.AddContextEventListeners();
		this.GetProgram();
	}

	static GetElement(id)
	{
		var element = document.getElementById(id);
		if(!element)
		{
			console.warn('No "' + id + '" element found');
		}
		return element;
	}
	static GetContext(canvasElement, contextId)
	{
		var context = canvasElement.getContext(contextId);
		if(!context)
		{
			console.warn('Could not get "' + id + '" context');
		}
		return context;
	}

	GetCanvasElement()
	{
		return this.canvasElement;
	}
	GetVertexBuffer()
	{
		if(!ShaderHelper.vertexBuffer)
		{
			ShaderHelper.vertexBuffer = this.context.createBuffer();
			this.context.bindBuffer(this.context.ARRAY_BUFFER,
				ShaderHelper.vertexBuffer);
			this.context.bufferData(this.context.ARRAY_BUFFER,
				new Float32Array(ShaderHelper.vertices), this.context.STATIC_DRAW);
			this.context.bindBuffer(this.context.ARRAY_BUFFER, null);
		}
		return ShaderHelper.vertexBuffer;
	}
	GetIndexBuffer()
	{
		if(!ShaderHelper.indexBuffer)
		{
			ShaderHelper.indexBuffer = this.context.createBuffer();
			this.context.bindBuffer(this.context.ELEMENT_ARRAY_BUFFER,
				ShaderHelper.indexBuffer);
			this.context.bufferData(this.context.ELEMENT_ARRAY_BUFFER,
				new Uint16Array(ShaderHelper.indices), this.context.STATIC_DRAW);
			this.context.bindBuffer(this.context.ELEMENT_ARRAY_BUFFER, null);
		}
		return ShaderHelper.indexBuffer;
	}
	GetVertexShader()
	{
		if(!this.vertexShader)
		{
			this.vertexShader = this.context.createShader(this.context.VERTEX_SHADER);
			this.context.shaderSource(this.vertexShader, this.vertexShaderCode);
			this.context.compileShader(this.vertexShader);
		}
		return this.vertexShader;
	}
	GetFragmentShader()
	{
		if(!this.fragmentShader)
		{
			this.fragmentShader = this.context.createShader(this.context.FRAGMENT_SHADER);
			this.context.shaderSource(this.fragmentShader, this.fragmentShaderCode);
			this.context.compileShader(this.fragmentShader);
		}
		return this.fragmentShader;
	}
	GetProgram()
	{
		if(!this.program)
		{
			this.program = this.context.createProgram();
			this.context.attachShader(this.program, this.GetVertexShader());
			this.context.attachShader(this.program, this.GetFragmentShader());

			this.context.linkProgram(this.program);
			this.context.useProgram(this.program);
			this.context.bindBuffer(this.context.ARRAY_BUFFER, this.GetVertexBuffer());
			this.context.bindBuffer(this.context.ELEMENT_ARRAY_BUFFER, this.GetIndexBuffer());

			var coord = this.context.getAttribLocation(this.program, "coordinates");
			this.context.vertexAttribPointer(coord, 3, this.context.FLOAT, false, 0, 0);
			this.context.enableVertexAttribArray(coord);

			this.context.viewport(0, 0, this.canvasElement.width, this.canvasElement.height);
		}
		return this.program;
	}
	DrawTriangles()
	{
		this.context.drawElements(this.context.TRIANGLES,
			ShaderHelper.indices.length, this.context.UNSIGNED_SHORT, 0);
	}
	DrawScene(now)
	{
		renderHelper.DrawTriangles();
		renderHelper.UpdateTime(now);
		renderHelper.DrawFrame();
	}
	DrawFrame()
	{
		this.SetUniform('width', this.canvasElement.width);
		this.SetUniform('height', this.canvasElement.height);
		this.SetUniform('time', 1000.0);
		this.requestId = requestAnimationFrame(this.DrawScene);
	}
	SetDrawTime(delta)
	{
		if(this.drawTimeElement)
		{
			this.drawTimeElement.textContent = delta.toFixed(3);
		}
	}
	Detach()
	{
		this.context.detachShader(this.program, this.fragmentShader);
		this.context.detachShader(this.program, this.vertexShader);
	}
	PrintLogs()
	{
		console.log('Shading language version: ' +
			this.context.getParameter(this.context.SHADING_LANGUAGE_VERSION));
		console.log('Fragment shader info log: ' +
			this.context.getShaderInfoLog(this.GetFragmentShader()));
		console.log('Vertex shader info log: ' +
			this.context.getShaderInfoLog(this.GetVertexShader()));
		console.log('Program info log: ' + this.context.getProgramInfoLog(this.program));
	}
	AddContextEventListeners()
	{
		this.canvasElement.addEventListener("webglcontextlost", this.OnContextLost, false);
		this.canvasElement.addEventListener("webglcontextrestored", this.OnContextRestored, false);
	}
	OnContextLost(event)
	{
		event.preventDefault();
		cancelAnimationFrame(this.requestId);
	}
	OnContextRestored(event)
	{
		//
	}
	SetUniform(id, value)
	{
		var handle = this.context.getUniformLocation(this.program, id);
		if(handle)
		{
			this.context.uniform1f(handle, value);
		}
	}
}
