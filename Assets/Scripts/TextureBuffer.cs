using UnityEngine;
using UnityEngine.Rendering;

public class TextureBuffer
{
    public RenderTexture[] TransmittanceArray { get; private set; }

    public RenderTexture[] IrradianceArray { get; private set; }

    public RenderTexture[] ScatteringArray { get; private set; }

    public RenderTexture[] OptionalSingleMieScatteringArray { get; private set; }

    public TextureBuffer(bool halfPrecision)
    {
        TransmittanceArray = NewTexture2DArray(
            CONSTANTS.TRANSMITTANCE_WIDTH, 
            CONSTANTS.TRANSMITTANCE_HEIGHT, 
            false);
    }

    public void Release()
    {
        ReleaseArray(TransmittanceArray);
    }

    public void Clear(ComputeShader compute)
    {
        ClearArray(compute, TransmittanceArray);
    }

    public static RenderTexture NewRenderTexture2D(int width, int height, bool halfPrecision)
    {
        RenderTextureFormat format = RenderTextureFormat.ARGBFloat;

        if (halfPrecision && SystemInfo.SupportsRenderTextureFormat(RenderTextureFormat.ARGBHalf))
            format = RenderTextureFormat.ARGBHalf;

        RenderTexture map = new RenderTexture(width, height, 0, format, RenderTextureReadWrite.Linear);
        map.filterMode = FilterMode.Bilinear;
        map.wrapMode = TextureWrapMode.Clamp;
        map.useMipMap = false;
        map.enableRandomWrite = true;
        map.Create();

        return map;
    }

    private RenderTexture[] NewTexture2DArray(int width, int height, bool halfPrecision)
    {
        RenderTexture[] arr = new RenderTexture[]
        {
            NewRenderTexture2D(width, height, halfPrecision),
            NewRenderTexture2D(width, height, halfPrecision)
        };
        return arr;
    }

    private void ReleaseArray(RenderTexture[] arr)
    {
        if (arr == null) return;

        for (int i = 0; i < arr.Length; i++)
        {
            ReleaseTexture(arr[i]);
            arr[i] = null;
        }
    }

    private void ReleaseTexture(RenderTexture tex)
    {
        if (tex == null) return;
        tex.Release();
        GameObject.DestroyImmediate(tex);
    }

    private void ClearArray(ComputeShader compute, RenderTexture[] arr)
    {
        if (arr == null) return;

        foreach (var tex in arr)
            ClearTexture(compute, tex);
    }

    private void ClearTexture(ComputeShader compute, RenderTexture tex)
    {
        if (tex == null) return;

        int NUM_THREADS = CONSTANTS.NUM_THREADS;

        if (tex.dimension == TextureDimension.Tex3D)
        {

        }
        else
        {
            int width = tex.width;
            int height = tex.height;

            int kernel = compute.FindKernel("ClearTex2D");
            compute.SetTexture(kernel, "targetWrite2D", tex);
            compute.Dispatch(kernel, width / NUM_THREADS, height / NUM_THREADS, 1);
        }
    }
}
