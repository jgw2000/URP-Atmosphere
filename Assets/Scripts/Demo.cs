using System;
using System.Collections.Generic;
using UnityEngine;

public class Demo : MonoBehaviour
{
    static readonly float kSunAngularRadius = 0.00935f / 2.0f;
    static readonly float kBottomRadius = 6360000.0f;
    static readonly float kLengthUnitInMeters = 1000.0f;

    public Light Sun;

    public bool UseConstantSolarSpectrum = false;

    public bool UseOzone = true;

    public bool UseCombinedTextures = true;

    public bool UseHalfPrecision = false;

    public bool DoWhiteBalance = false;

    public LUMINANCE UseLuminance = LUMINANCE.NONE;

    public float Exposure = 10.0f;

    public ComputeShader m_compute;

    public Material m_material;

    private Model m_model;

    void Awake()
    {
        int kLambdaMin = 360;
        int kLambdaMax = 830;

        // 太阳辐照数据，摘自ASTM G-173 标准光谱（http://rredc.nrel.gov/solar/spectra/am1.5/ASTMG173/ASTMG173.html）
        // 用于精确描述太阳光谱在地表附近（标准大气压下）的分布，数组中的每个值代表一个波段（10nm）内的太阳辐照度，
        // 单位是 W.m^-2
        double[] kSolarIrradiance = new double[]
        {
            1.11776, 1.14259, 1.01249, 1.14716, 1.72765, 1.73054, 1.6887, 1.61253,
            1.91198, 2.03474, 2.02042, 2.02212, 1.93377, 1.95809, 1.91686, 1.8298,
            1.8685, 1.8931, 1.85149, 1.8504, 1.8341, 1.8345, 1.8147, 1.78158, 1.7533,
            1.6965, 1.68194, 1.64654, 1.6048, 1.52143, 1.55622, 1.5113, 1.474, 1.4482,
            1.41018, 1.36775, 1.34188, 1.31429, 1.28303, 1.26758, 1.2367, 1.2082,
            1.18737, 1.14683, 1.12362, 1.1058, 1.07124, 1.04992
        };

        // 臭氧截面数据，参考 http://www.iup.uni-bremen.de/gruppen/molspec/databases/referencespectra/o3spectra2011/index.html
        // 数组中的每个值代表一个波段（10nm）内每个臭氧分子的光学截面，单位是 m^2
        double[] kOzoneCrossSection = new double[]
        {
            1.18e-27, 2.182e-28, 2.818e-28, 6.636e-28, 1.527e-27, 2.763e-27, 5.52e-27,
            8.451e-27, 1.582e-26, 2.316e-26, 3.669e-26, 4.924e-26, 7.752e-26, 9.016e-26,
            1.48e-25, 1.602e-25, 2.139e-25, 2.755e-25, 3.091e-25, 3.5e-25, 4.266e-25,
            4.672e-25, 4.398e-25, 4.701e-25, 5.019e-25, 4.305e-25, 3.74e-25, 3.215e-25,
            2.662e-25, 2.238e-25, 1.852e-25, 1.473e-25, 1.209e-25, 9.423e-26, 7.455e-26,
            6.566e-26, 5.105e-26, 4.15e-26, 4.228e-26, 3.237e-26, 2.451e-26, 2.801e-26,
            2.534e-26, 1.624e-26, 1.465e-26, 2.078e-26, 1.383e-26, 7.105e-27
         };

        // 参考  https://en.wikipedia.org/wiki/Dobson_unit，表示 1 DU 的臭氧分子的数量，单位是 m^-2
        double kDobsonUnit = 2.687e20;

        // 最大臭氧分子密度
        double kMaxOzoneNumberDensity = 300.0 * kDobsonUnit / 15000.0;

        double kConstantSolarIrradiance = 1.5;
        double kTopRadius = 6420000.0;
        double kRayleigh = 1.24062e-6;
        double kRayleighScaleHeight = 8000.0;
        double kMieScaleHeight = 1200.0;
        double kMieAngstromAlpha = 0.0;
        double kMieAngstromBeta = 5.328e-3;
        double kMieSingleScatteringAlbedo = 0.9;
        double kMiePhaseFunctionG = 0.8;
        double kGroundAlbedo = 0.1;
        double max_sun_zenith_angle = (UseHalfPrecision ? 102.0 : 120.0) / 180.0 * Mathf.PI;

        DensityProfileLayer rayleigh_layer = new DensityProfileLayer("rayleigh", 0.0, 1.0, -1.0 / kRayleighScaleHeight, 0.0, 0.0);
        DensityProfileLayer mie_layer = new DensityProfileLayer("mie", 0.0, 1.0, -1.0 / kMieScaleHeight, 0.0, 0.0);

        // 臭氧层 10-25 km 之间从 0 线性增加到 1，25-40 km 之间从 1 线性降低到 0
        List<DensityProfileLayer> ozone_density = new List<DensityProfileLayer>();
        ozone_density.Add(new DensityProfileLayer("absorption0", 25000.0, 0.0, 0.0, 1.0 / 15000.0, -2.0 / 3.0));
        ozone_density.Add(new DensityProfileLayer("absorption1", 0.0, 0.0, 0.0, -1.0 / 15000.0, 8.0 / 3.0));

        List<double> wavelengths = new List<double>();
        List<double> solar_irradiance = new List<double>();
        List<double> rayleigh_scattering = new List<double>();
        List<double> mie_scattering = new List<double>();
        List<double> mie_extinction = new List<double>();
        List<double> absorption_extinction = new List<double>();
        List<double> ground_albedo = new List<double>();

        for (int l = kLambdaMin; l <= kLambdaMax; l += 10)
        {
            double lambda = l * 1e-3;   // 微米
            double mie = kMieAngstromBeta / kMieScaleHeight * Math.Pow(lambda, -kMieAngstromAlpha);

            wavelengths.Add(l);

            if (UseConstantSolarSpectrum)
                solar_irradiance.Add(kConstantSolarIrradiance);
            else
                solar_irradiance.Add(kSolarIrradiance[(l - kLambdaMin) / 10]);

            rayleigh_scattering.Add(kRayleigh * Math.Pow(lambda, -4));
            mie_scattering.Add(mie * kMieSingleScatteringAlbedo);
            mie_extinction.Add(mie);
            absorption_extinction.Add(UseOzone ? kMaxOzoneNumberDensity * kOzoneCrossSection[(l - kLambdaMin) / 10] : 0.0);
            ground_albedo.Add(kGroundAlbedo);
        }

        m_model = new Model();

        m_model.HalfPrecision = UseHalfPrecision;
        m_model.CombineScatteringTextures = UseCombinedTextures;
        m_model.UseLuminance = UseLuminance;
        m_model.Wavelengths = wavelengths;
        m_model.SolarIrradiance = solar_irradiance;
        m_model.BottomRadius = kBottomRadius;
        m_model.TopRadius = kTopRadius;
        m_model.RayleighDensity = rayleigh_layer;
        m_model.RayleighScattering = rayleigh_scattering;
        m_model.MieDensity = mie_layer;
        m_model.MieScattering = mie_scattering;
        m_model.MieExtinction = mie_extinction;
        m_model.MiePhaseFunctionG = kMiePhaseFunctionG;
        m_model.AbsorptionDensity = ozone_density;
        m_model.AbsorptionExtinction = absorption_extinction;
        m_model.GroundAlbedo = ground_albedo;
        m_model.MaxSunZenithAngle = max_sun_zenith_angle;
        m_model.LengthUnitInMeters = kLengthUnitInMeters;

        int numScatteringOrders = 4;
        m_model.Init(m_compute, numScatteringOrders);
    }

    private void OnDestroy()
    {
        m_model?.Release();
    }
}
