
class RenderHelper
{
	#previousFrameTime = 0;
	#shaderHelper = null;

	constructor(shaderHelper)
	{
		this.SetShaderHelper(shaderHelper);
	}

	UpdateTime(now)
	{
		this.shaderHelper.SetDrawTime(now - this.previousFrameTime);
		this.previousFrameTime = now;
	}
	GetShaderHelper()
	{
		if(!this.shaderHelper)
		{
			console.log('No shader helper added');
		}
		return this.shaderHelper;
	}
	SetShaderHelper(shaderHelper)
	{
		if(this.shaderHelper)
		{
			this.shaderHelper.Detach();
		}
		this.shaderHelper = shaderHelper;
	}
	GetProgram()
	{
		this.shaderHelper.GetProgram(this.canvasElement);
	}
	PrintLogs()
	{
		this.shaderHelper.PrintLogs();
	}
	DrawTriangles()
	{
		this.shaderHelper.DrawTriangles();
	}
	DrawFrame()
	{
		this.shaderHelper.DrawFrame();
	}
}
