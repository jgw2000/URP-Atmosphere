using UnityEngine;
using UnityEngine.Rendering;

public class TextureBuffer
{
    public RenderTexture DeltaIrradianceTexture {  get; private set; }

    public RenderTexture DeltaRayleighScatteringTexture { get; private set; }

    public RenderTexture DeltaMieScatteringTexture { get; private set; }

    public RenderTexture DeltaScatteringDensityTexture { get; private set; }

    public RenderTexture DeltaMultipleScatteringTexture { get; private set; }

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

        IrradianceArray = NewTexture2DArray(
            CONSTANTS.IRRADIANCE_WIDTH,
            CONSTANTS.IRRADIANCE_HEIGHT,
            false);

        ScatteringArray = NewTexture3DArray(
            CONSTANTS.SCATTERING_WIDTH,
            CONSTANTS.SCATTERING_HEIGHT,
            CONSTANTS.SCATTERING_DEPTH,
            halfPrecision);

        OptionalSingleMieScatteringArray = NewTexture3DArray(
            CONSTANTS.SCATTERING_WIDTH,
            CONSTANTS.SCATTERING_HEIGHT,
            CONSTANTS.SCATTERING_DEPTH,
            halfPrecision);

        DeltaIrradianceTexture = NewRenderTexture2D(
            CONSTANTS.IRRADIANCE_WIDTH,
            CONSTANTS.IRRADIANCE_HEIGHT,
            false);

        DeltaRayleighScatteringTexture = NewRenderTexture3D(
            CONSTANTS.SCATTERING_WIDTH,
            CONSTANTS.SCATTERING_HEIGHT,
            CONSTANTS.SCATTERING_DEPTH,
            halfPrecision);

        DeltaMieScatteringTexture = NewRenderTexture3D(
            CONSTANTS.SCATTERING_WIDTH,
            CONSTANTS.SCATTERING_HEIGHT,
            CONSTANTS.SCATTERING_DEPTH,
            halfPrecision);

        DeltaScatteringDensityTexture = NewRenderTexture3D(
            CONSTANTS.SCATTERING_WIDTH,
            CONSTANTS.SCATTERING_HEIGHT,
            CONSTANTS.SCATTERING_DEPTH,
            halfPrecision);

        // delta_multiple_scattering_texture is only needed to compute scattering
        // order 3 or more, while delta_rayleigh_scattering_texture and
        // delta_mie_scattering_texture are only needed to compute double scattering.
        // Therefore, to save memory, we can store delta_rayleigh_scattering_texture
        // and delta_multiple_scattering_texture in the same GPU texture.
        DeltaMultipleScatteringTexture = DeltaRayleighScatteringTexture;
    }

    public void Release()
    {
        ReleaseTexture(DeltaIrradianceTexture);
        ReleaseTexture(DeltaRayleighScatteringTexture);
        ReleaseTexture(DeltaMieScatteringTexture);
        ReleaseTexture(DeltaScatteringDensityTexture);
        ReleaseArray(TransmittanceArray);
        ReleaseArray(IrradianceArray);
        ReleaseArray(ScatteringArray);
        ReleaseArray(OptionalSingleMieScatteringArray);
    }

    public void Clear(ComputeShader compute)
    {
        ClearTexture(compute, DeltaIrradianceTexture);
        ClearTexture(compute, DeltaRayleighScatteringTexture);
        ClearTexture(compute, DeltaMieScatteringTexture);
        ClearTexture(compute, DeltaScatteringDensityTexture);
        ClearArray(compute, TransmittanceArray);
        ClearArray(compute, IrradianceArray);
        ClearArray(compute, ScatteringArray);
        ClearArray(compute, OptionalSingleMieScatteringArray);
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

    public static RenderTexture NewRenderTexture3D(int width, int height, int depth, bool halfPrecision)
    {
        RenderTextureFormat format = RenderTextureFormat.ARGBFloat;

        if (halfPrecision && SystemInfo.SupportsRenderTextureFormat(RenderTextureFormat.ARGBHalf))
            format = RenderTextureFormat.ARGBHalf;

        RenderTexture map = new RenderTexture(width, height, 0, format, RenderTextureReadWrite.Linear);
        map.volumeDepth = depth;
        map.dimension = TextureDimension.Tex3D;
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

    private RenderTexture[] NewTexture3DArray(int width, int height, int depth, bool halfPrecision)
    {
        RenderTexture[] arr = new RenderTexture[]
        {
            NewRenderTexture3D(width, height, depth, halfPrecision),
            NewRenderTexture3D(width, height, depth, halfPrecision)
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
            int width = tex.width;
            int height = tex.height;
            int depth = tex.volumeDepth;

            int kernel = compute.FindKernel("ClearTex3D");
            compute.SetTexture(kernel, "targetWrite3D", tex);
            compute.Dispatch(kernel, width / NUM_THREADS, height / NUM_THREADS, depth / NUM_THREADS);
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
