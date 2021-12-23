
const renderHelper = new RenderHelper(shader_basic);
renderHelper.DrawFrame();
renderHelper.PrintLogs();
renderHelper.SetShaderHelper(shader_colors);

// renderHelper.GetShaderHelper().SetUniform('width',
// 	renderHelper.GetShaderHelper().GetCanvasElement().width);
// renderHelper.GetShaderHelper().SetUniform('height',
// 	renderHelper.GetShaderHelper().GetCanvasElement().height);
// renderHelper.GetShaderHelper().SetUniform('time', 1000.0);

renderHelper.PrintLogs();
