using UnityEngine;
using UnityEngine.Rendering.Universal;

public class AtmosphereFeature : ScriptableRendererFeature
{
    public Material atmosphereMaterial;
    public ComputeShader computeShader;
    public RenderTexture transmittanceLUT;

    private AtmospherePass atmospherePass;

    public override void Create()
    {
        atmospherePass = new AtmospherePass(atmosphereMaterial);
    }

    public override void SetupRenderPasses(ScriptableRenderer renderer, in RenderingData renderingData)
    {
        atmospherePass.SetTarget(renderer.cameraColorTargetHandle);
    }

    public override void AddRenderPasses(ScriptableRenderer renderer, ref RenderingData renderingData)
    {
        if (atmospherePass != null)
        {
            renderer.EnqueuePass(atmospherePass);
        }
    }

    public void BakeTransmittanceLUT()
    {
        computeShader.SetTexture(0, "_Result", transmittanceLUT);
        computeShader.SetInt("_TextureSize", transmittanceLUT.width);

        computeShader.GetKernelThreadGroupSizes(0, out uint x, out uint y, out _);
        int numGroupsX = Mathf.CeilToInt((float)transmittanceLUT.width / x);
        int numGroupsY = Mathf.CeilToInt((float)transmittanceLUT.height / y);

        computeShader.Dispatch(0, numGroupsY, numGroupsY, 1);
    }
}
