using System;
using UnityEngine;
using UnityEngine.Rendering;
using UnityEngine.Rendering.Universal;

public class AtmospherePass : ScriptableRenderPass
{
    private ProfilingSampler m_ProfilingSampler = new ProfilingSampler("Atmosphere");
    private Material m_AtmosphereMaterial;
    private RTHandle m_CameraColorTarget;

    public AtmospherePass(Material material)
    {
        this.m_AtmosphereMaterial = material;
        renderPassEvent = RenderPassEvent.BeforeRenderingOpaques;
    }

    public void SetTarget(RTHandle colorHandle)
    {
        m_CameraColorTarget = colorHandle;
    }

    [Obsolete]
    public override void Execute(ScriptableRenderContext context, ref RenderingData renderingData)
    {
        var cameraData = renderingData.cameraData;
        if (cameraData.camera.cameraType != CameraType.SceneView && cameraData.camera.cameraType != CameraType.Game)
            return;

        if (m_AtmosphereMaterial == null)
        {
            return;
        }

        CommandBuffer cmd = CommandBufferPool.Get("Atmosphere Pass");
        using (new ProfilingScope(cmd, m_ProfilingSampler))
        {
            Blitter.BlitCameraTexture(cmd, m_CameraColorTarget, m_CameraColorTarget, m_AtmosphereMaterial, 0);
        }
        context.ExecuteCommandBuffer(cmd);
        cmd.Clear();

        CommandBufferPool.Release(cmd);
    }
}
