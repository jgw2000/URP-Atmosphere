using System;
using System.Collections.Generic;
using UnityEngine;

public class Model
{
    private const int READ = 0;
    private const int WRITE = 1;

    private const double kLambdaR = 680.0;
    private const double kLambdaG = 550.0;
    private const double kLambdaB = 440.0;

    private const int kLambdaMin = 360;
    private const int kLambdaMax = 830;

    public IList<double> Wavelengths { get; set; }

    public IList<double> SolarIrradiance { get; set; }

    public double SunAngularRadius { get; set; }

    public double BottomRadius { get; set; }

    public double TopRadius { get; set; }

    public DensityProfileLayer RayleighDensity { get; set; }

    public IList<double> RayleighScattering { get; set; }

    public DensityProfileLayer MieDensity { get; set; }

    public IList<double> MieScattering { get; set; }

    public IList<double> MieExtinction { get; set; }

    public double MiePhaseFunctionG { get; set; }

    public IList<DensityProfileLayer> AbsorptionDensity { get; set; }

    public IList<double> AbsorptionExtinction { get; set; }

    public IList<double> GroundAlbedo {  get; set; }

    /// <summary>
    /// The maximum Sun zenith angle for which atmospheric scattering must be
    /// precomputed, in radians (for maximum precision, use the smallest Sun
    /// zenith angle yielding negligible sky light radiance values. For instance,
    /// for the Earth case, 102 degrees is a good choice for most cases (120
    /// degrees is necessary for very high exposure values).
    /// </summary>
    public double MaxSunZenithAngle { get; set; }

    /// <summary>
    /// The length unit used in your shaders and meshes. This is the length unit
    /// which must be used when calling the atmosphere model shader functions.
    /// </summary>
    public double LengthUnitInMeters { get; set; }

    public int NumPrecomputedWavelengths { get { return UseLuminance == LUMINANCE.PRECOMPUTED ? 15 : 3; } }

    public bool CombineScatteringTextures { get; set; }

    public LUMINANCE UseLuminance { get; set; }

    public bool HalfPrecision {  get; set; }

    public RenderTexture TransmittanceTexture { get; private set; }

    public RenderTexture ScatteringTexture { get; private set; }

    public RenderTexture IrradianceTexture { get; private set; }

    public RenderTexture OptionalSingleMieScatteringTexture { get; private set; }

    public Model()
    {

    }

    public void Init(ComputeShader compute, int num_scattering_orders)
    {
        TextureBuffer buffer = new TextureBuffer(HalfPrecision);
        buffer.Clear(compute);

        if (NumPrecomputedWavelengths <= 3)
        {
            Precompute(compute, buffer, null, null, false, num_scattering_orders);
        }

        TransmittanceTexture = buffer.TransmittanceArray[READ];
        buffer.TransmittanceArray[READ] = null;

        ScatteringTexture = buffer.ScatteringArray[READ];
        buffer.ScatteringArray[READ] = null;

        buffer.Release();
    }

    public void Release()
    {
        ReleaseTexture(TransmittanceTexture);
    }

    void Precompute(
        ComputeShader compute,
        TextureBuffer buffer,
        double[] lambdas,
        double[] luminance_from_radiance,
        bool blend,
        int num_scattering_orders)
    {
        int BLEND = blend ? 1 : 0;
        int NUM_THREADS = CONSTANTS.NUM_THREADS;

        BindToCompute(compute, lambdas, luminance_from_radiance);

        int compute_transmittance = compute.FindKernel("ComputeTransmittance");
        int compute_single_scattering = compute.FindKernel("ComputeSingleScattering");

        // Compute the transmittance, and store it in transmittance_texture
        compute.SetTexture(compute_transmittance, "transmittanceWrite", buffer.TransmittanceArray[WRITE]);
        compute.SetVector("blend", new Vector4(0, 0, 0, 0));
        compute.Dispatch(compute_transmittance, CONSTANTS.TRANSMITTANCE_WIDTH / NUM_THREADS, CONSTANTS.TRANSMITTANCE_HEIGHT / NUM_THREADS, 1);
        Swap(buffer.TransmittanceArray);

        // Compute the rayleigh and mie single scattering
        compute.SetTexture(compute_single_scattering, "scatteringWrite", buffer.ScatteringArray[WRITE]);
        compute.SetTexture(compute_single_scattering, "_transmittance_texture", buffer.TransmittanceArray[READ]);
        compute.SetVector("blend", new Vector4(0, 0, BLEND, BLEND));
        compute.Dispatch(compute_single_scattering, CONSTANTS.SCATTERING_WIDTH / NUM_THREADS, CONSTANTS.SCATTERING_HEIGHT / NUM_THREADS, CONSTANTS.SCATTERING_DEPTH /  NUM_THREADS);
        Swap(buffer.ScatteringArray);
    }

    public void BindToMaterial(Material mat)
    {
        mat.SetTexture("_transmittance_texture", TransmittanceTexture);
        mat.SetTexture("_scattering_texture", ScatteringTexture);

        mat.SetInt("TRANSMITTANCE_TEXTURE_WIDTH", CONSTANTS.TRANSMITTANCE_WIDTH);
        mat.SetInt("TRANSMITTANCE_TEXTURE_HEIGHT", CONSTANTS.TRANSMITTANCE_HEIGHT);
        mat.SetInt("SCATTERING_TEXTURE_R_SIZE", CONSTANTS.SCATTERING_R);
        mat.SetInt("SCATTERING_TEXTURE_MU_SIZE", CONSTANTS.SCATTERING_MU);
        mat.SetInt("SCATTERING_TEXTURE_MU_S_SIZE", CONSTANTS.SCATTERING_MU_S);
        mat.SetInt("SCATTERING_TEXTURE_NU_SIZE", CONSTANTS.SCATTERING_NU);
        mat.SetInt("SCATTERING_TEXTURE_WIDTH", CONSTANTS.SCATTERING_WIDTH);
        mat.SetInt("SCATTERING_TEXTURE_HEIGHT", CONSTANTS.SCATTERING_HEIGHT);
        mat.SetInt("SCATTERING_TEXTURE_DEPTH", CONSTANTS.SCATTERING_DEPTH);

        mat.SetFloat("sun_angular_radius", (float)SunAngularRadius);
        mat.SetFloat("bottom_radius", (float)(BottomRadius / LengthUnitInMeters));
        mat.SetFloat("top_radius", (float)(TopRadius / LengthUnitInMeters));
        mat.SetFloat("mie_phase_function_g", (float)MiePhaseFunctionG);
        mat.SetFloat("mu_s_min", (float)Math.Cos(MaxSunZenithAngle));

        BindDensityLayer(mat, RayleighDensity);
        BindDensityLayer(mat, MieDensity);

        double[] lambdas = new double[] { kLambdaR, kLambdaG, kLambdaB };

        Vector3 solarIrradiance = ToVector(Wavelengths, SolarIrradiance, lambdas, 1.0);
        mat.SetVector("solar_irradiance", solarIrradiance);

        Vector3 rayleighScattering = ToVector(Wavelengths, RayleighScattering, lambdas, LengthUnitInMeters);
        mat.SetVector("rayleigh_scattering", rayleighScattering);

        Vector3 mieScattering = ToVector(Wavelengths, MieScattering, lambdas, LengthUnitInMeters);
        mat.SetVector("mie_scattering", mieScattering);
    }

    /// <summary>
    /// Bind to a compute shader for precomutation of textures.
    /// </summary>
    private void BindToCompute(ComputeShader compute, double[] lambdas, double[] luminance_from_radiance)
    {
        if (lambdas == null)
            lambdas = new double[] { kLambdaR, kLambdaG, kLambdaB };

        if (luminance_from_radiance == null)
            luminance_from_radiance = new double[] { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };

        compute.SetInt("TRANSMITTANCE_TEXTURE_WIDTH", CONSTANTS.TRANSMITTANCE_WIDTH);
        compute.SetInt("TRANSMITTANCE_TEXTURE_HEIGHT", CONSTANTS.TRANSMITTANCE_HEIGHT);
        compute.SetInt("SCATTERING_TEXTURE_R_SIZE", CONSTANTS.SCATTERING_R);
        compute.SetInt("SCATTERING_TEXTURE_MU_SIZE", CONSTANTS.SCATTERING_MU);
        compute.SetInt("SCATTERING_TEXTURE_MU_S_SIZE", CONSTANTS.SCATTERING_MU_S);
        compute.SetInt("SCATTERING_TEXTURE_NU_SIZE", CONSTANTS.SCATTERING_NU);
        compute.SetInt("SCATTERING_TEXTURE_WIDTH", CONSTANTS.SCATTERING_WIDTH);
        compute.SetInt("SCATTERING_TEXTURE_HEIGHT", CONSTANTS.SCATTERING_HEIGHT);
        compute.SetInt("SCATTERING_TEXTURE_DEPTH", CONSTANTS.SCATTERING_DEPTH);

        Vector3 solarIrradiance = ToVector(Wavelengths, SolarIrradiance, lambdas, 1.0);
        compute.SetVector("solar_irradiance", solarIrradiance);

        Vector3 rayleighScattering = ToVector(Wavelengths, RayleighScattering, lambdas, LengthUnitInMeters);
        BindDensityLayer(compute, RayleighDensity);
        compute.SetVector("rayleigh_scattering", rayleighScattering);

        Vector3 mieScattering = ToVector(Wavelengths, MieScattering, lambdas, LengthUnitInMeters);
        Vector3 mieExtinction = ToVector(Wavelengths, MieExtinction, lambdas, LengthUnitInMeters);
        BindDensityLayer(compute, MieDensity);
        compute.SetVector("mie_scattering", mieScattering);
        compute.SetVector("mie_extinction", mieExtinction);

        Vector3 absorptionExtinction = ToVector(Wavelengths, AbsorptionExtinction, lambdas, LengthUnitInMeters);
        BindDensityLayer(compute, AbsorptionDensity[0]);
        BindDensityLayer(compute, AbsorptionDensity[1]);
        compute.SetVector("absorption_extinction", absorptionExtinction);

        compute.SetFloat("sun_angular_radius", (float)SunAngularRadius);
        compute.SetFloat("bottom_radius", (float)(BottomRadius / LengthUnitInMeters));
        compute.SetFloat("top_radius", (float)(TopRadius / LengthUnitInMeters));
        compute.SetFloat("mie_phase_function_g", (float)MiePhaseFunctionG);
        compute.SetFloat("mu_s_min", (float)Math.Cos(MaxSunZenithAngle));
    }

    private void BindDensityLayer(Material mat, DensityProfileLayer layer)
    {
        mat.SetFloat(layer.Name + "_width", (float)(layer.Width / LengthUnitInMeters));
        mat.SetFloat(layer.Name + "_exp_term", (float)layer.ExpTerm);
        mat.SetFloat(layer.Name + "_exp_scale", (float)(layer.ExpScale * LengthUnitInMeters));
        mat.SetFloat(layer.Name + "_linear_term", (float)(layer.LinearTerm * LengthUnitInMeters));
        mat.SetFloat(layer.Name + "_constant_term", (float)layer.ConstantTerm);
    }

    private void BindDensityLayer(ComputeShader compute, DensityProfileLayer layer)
    {
        compute.SetFloat(layer.Name + "_width", (float)(layer.Width / LengthUnitInMeters));
        compute.SetFloat(layer.Name + "_exp_term", (float)layer.ExpTerm);
        compute.SetFloat(layer.Name + "_exp_scale", (float)(layer.ExpScale * LengthUnitInMeters));
        compute.SetFloat(layer.Name + "_linear_term", (float)(layer.LinearTerm * LengthUnitInMeters));
        compute.SetFloat(layer.Name + "_constant_term", (float)layer.ConstantTerm);
    }

    private Vector3 ToVector(IList<double> wavelengths, IList<double> v, IList<double> lambdas, double scale)
    {
        double r = Interpolate(wavelengths, v, lambdas[0]) * scale;
        double g = Interpolate(wavelengths, v, lambdas[1]) * scale;
        double b = Interpolate(wavelengths, v, lambdas[2]) * scale;

        return new Vector3((float)r, (float)g, (float)b);
    }

    private static double Interpolate(IList<double> wavelengths, IList<double> wavelength_function, double wavelength)
    {
        if (wavelength < wavelengths[0]) return wavelength_function[0];

        for (int i = 0; i < wavelengths.Count - 1; ++i)
        {
            if (wavelength < wavelengths[i + 1])
            {
                double u = (wavelength - wavelengths[i]) / (wavelengths[i + 1] - wavelengths[i]);
                return wavelength_function[i] * (1.0 - u) + wavelength_function[i + 1] * u;
            }
        }

        return wavelength_function[wavelength_function.Count - 1];
    }

    private void Swap(RenderTexture[] arr)
    {
        RenderTexture tmp = arr[READ];
        arr[READ] = arr[WRITE];
        arr[WRITE] = tmp;
    }

    private void ReleaseTexture(RenderTexture tex)
    {
        if (tex == null) return;
        tex.Release();
        GameObject.DestroyImmediate(tex);
    }
}
