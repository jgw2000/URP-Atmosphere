#pragma once

#include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/Lighting.hlsl"
#include "Common.hlsl"

static const int _InScatteringPoints = 20;

CBUFFER_START(UnityPerMaterial)
float _MieG;
CBUFFER_END

TEXTURE2D(_TransmittanceLUT);
SAMPLER(sampler_TransmittanceLUT);

float3 GetTransmittanceToTopAtmosphereBoundary(float r, float mu)
{
    float2 uv = GetTransmittanceTextureUvFromRMu(r, mu);
    return SAMPLE_TEXTURE2D(_TransmittanceLUT, sampler_TransmittanceLUT, uv).xyz;
}

float3 GetTransmittance(float r, float mu, float d, bool ray_r_mu_intersects_ground)
{
    float r_d = ClampRadius(sqrt(d * d + 2.0 * r * mu * d + r * r));
    float mu_d = ClampCosine((r * mu + d) / r_d);
    
    if (ray_r_mu_intersects_ground)
    {
        return min(
            GetTransmittanceToTopAtmosphereBoundary(r_d, -mu_d) / GetTransmittanceToTopAtmosphereBoundary(r, -mu),
            1.0
        );
    }
    else
    {
        return min(
            GetTransmittanceToTopAtmosphereBoundary(r, mu) / GetTransmittanceToTopAtmosphereBoundary(r_d, mu_d),
            1.0
        );
    }
}

float3 GetTransmittanceToSun(float r, float mu_s)
{
    float sin_theta_h = kBottomRadius / r;
    float cos_theta_h = -sqrt(max(1.0 - sin_theta_h * sin_theta_h, 0.0));
    return GetTransmittanceToTopAtmosphereBoundary(r, mu_s) * smoothstep(-sin_theta_h * kSunAngularRadius, sin_theta_h * kSunAngularRadius, mu_s - cos_theta_h);
}

float3 ComputeSingleScatteringIntegrand(float r, float mu, float mu_s, float nu, float d, bool ray_r_mu_intersects_ground)
{
    float r_d = ClampRadius(sqrt(d * d + 2.0 * r * mu * d + r * r));
    float mu_s_d = ClampCosine((r * mu_s + d * nu) / r_d);
    float3 transmittance = GetTransmittance(r, mu, d, ray_r_mu_intersects_ground) * GetTransmittanceToSun(r_d, mu_s_d);
    float3 rayleigh = transmittance * GetLayerDensity(kRayleighScaleHeight, r_d - kBottomRadius);
    return rayleigh;
}

float DistanceToNearestAtmosphereBoundary(float r, float mu, bool ray_r_mu_intersects_ground)
{
    if (ray_r_mu_intersects_ground)
    {
        return DistanceToBottomAtmosphereBoundary(r, mu);
    }
    else
    {
        return DistanceToTopAtmosphereBoundary(r, mu);
    }
}

float3 ComputeSingleScattering(float r, float mu, float mu_s, float nu, bool ray_r_mu_intersects_ground)
{
    const int SAMPLE_COUNT = 50;
    float dx = DistanceToNearestAtmosphereBoundary(r, mu, ray_r_mu_intersects_ground);
    float3 rayleigh_sum = 0.0;
    
    for (int i = 0; i <= SAMPLE_COUNT; ++i)
    {
        float d_i = i * dx;
        float3 rayleigh_i = ComputeSingleScatteringIntegrand(r, mu, mu_s, nu, d_i, ray_r_mu_intersects_ground);
        float weight_i = (i == 0 || i == SAMPLE_COUNT) ? 0.5 : 1.0;
        rayleigh_sum += rayleigh_i * weight_i;
    }
    
    Light sun = GetMainLight();
    float3 rayleigh = rayleigh_sum * dx * sun.color * kRayleighScattering;
    return rayleigh;
}

float3 GetSkyRadiance(float3 p, float3 view_ray)
{
    float r = length(p);
    float rmu = dot(p, view_ray);
    
    Light sun = GetMainLight();
    
    float mu = rmu / r;
    float mu_s = dot(p, sun.direction) / r;
    float nu = dot(view_ray, sun.direction);
    bool ray_r_mu_intersects_ground = RayIntersectsGround(r, mu);
    
    return ComputeSingleScattering(r, mu, mu_s, nu, ray_r_mu_intersects_ground);
}

float3 OpticalDepthBaked(float3 sphereCenter, float3 rayOrigin, float3 sunDir)
{
    float rayLen = length(rayOrigin - sphereCenter);
    float h = rayLen - Rg;
    float uvY = saturate(h / (Rt - Rg));
    
    float3 rayDir = (rayOrigin - sphereCenter) / rayLen;
    float uvX = 1 - (dot(rayDir, sunDir) * 0.5 + 0.5);
    
    return SAMPLE_TEXTURE2D(_TransmittanceLUT, sampler_TransmittanceLUT, float2(uvX, uvY)).xyz;
}

float3 CalculateAtmosphere(float3 rayOrigin, float3 rayDir)
{
    float3 sphereCenter = float3(0, 0, 0);
    float distance = RaySphere(sphereCenter, Rt, rayOrigin, rayDir);
    float stepSize = distance / _InScatteringPoints;
    float3 inScatterPoint = rayOrigin + rayDir * stepSize * 0.5;
    
    Light sun = GetMainLight();
    
    float3 totalRay = 0;
    float3 totalMie = 0;
    float3 opticalDepth = 0;
    
    [unroll(MAX_LOOP_ITERATIONS)]
    for (int i = 0; i < _InScatteringPoints; i++)
    {
        // Particle density at sample position
        float3 density = DensityAtPoint(sphereCenter, inScatterPoint) * stepSize;
        
        // Accumulate optical depth
        opticalDepth += density;
        
        // Light ray optical depth
        float3 lightOpticalDepth = OpticalDepthBaked(sphereCenter, inScatterPoint, sun.direction);
        
        // Attenuation calculation
        float3 attenuation = exp(
            -betaR * (opticalDepth.x + lightOpticalDepth.x)
            -betaM * (opticalDepth.y + lightOpticalDepth.y)
        );
        
        // Accumulate scattered light
        totalRay += density.x * attenuation;
        totalMie += density.y * attenuation;
        
        inScatterPoint += rayDir * stepSize;
    }
    
    float mu = dot(rayDir, normalize(sun.direction));
    float phaseRay = 3.0 / (16.0 * PI) * (1 + mu * mu);
    float phaseMie = 3.0 / (8.0 * PI) * ((1 - _MieG * _MieG) * (1 + mu * mu)) / (pow(abs(1 + _MieG * _MieG - 2 * _MieG * mu), 1.5) * (2.0 + _MieG * _MieG));
    
    // Calculate final scattering factors
    float3 rayleigh = phaseRay * betaR * totalRay;
    float3 mie = phaseMie * betaM * totalMie;
    
    return sun.color * (rayleigh + mie);
}