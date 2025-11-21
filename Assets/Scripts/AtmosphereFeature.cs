using UnityEngine;
using UnityEngine.Rendering.Universal;

public class AtmosphereFeature : ScriptableRendererFeature
{
    public Material atmosphereMaterial;

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
}
