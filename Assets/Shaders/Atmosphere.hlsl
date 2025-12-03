#pragma once

#include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/Lighting.hlsl"
#include "Common.hlsl"

CBUFFER_START(UnityPerMaterial)
float _MieG;
CBUFFER_END

TEXTURE2D(_TransmittanceLUT);
SAMPLER(sampler_TransmittanceLUT);

TEXTURE3D(_ScatteringLUT);
SAMPLER(sampler_ScatteringLUT);

float3 GetTransmittanceToTopAtmosphereBoundary(float r, float mu)
{
    float2 uv = GetTransmittanceTextureUvFromRMu(r, mu);
    return SAMPLE_TEXTURE2D(_TransmittanceLUT, sampler_TransmittanceLUT, uv).rgb;
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

void ComputeSingleScatteringIntegrand(float r, float mu, float mu_s, float nu, float d, bool ray_r_mu_intersects_ground, out float3 rayleigh, out float3 mie)
{
    float r_d = ClampRadius(sqrt(d * d + 2.0 * r * mu * d + r * r));
    float mu_s_d = ClampCosine((r * mu_s + d * nu) / r_d);
    float3 transmittance = GetTransmittance(r, mu, d, ray_r_mu_intersects_ground) * GetTransmittanceToSun(r_d, mu_s_d);
    rayleigh = transmittance * GetLayerDensity(kRayleighScaleHeight, r_d - kBottomRadius);
    mie = transmittance * GetLayerDensity(kMieScaleHeight, r_d - kBottomRadius);
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

void ComputeSingleScattering(float r, float mu, float mu_s, float nu, bool ray_r_mu_intersects_ground, out float3 rayleigh, out float3 mie)
{
    const int SAMPLE_COUNT = 50;
    float dx = DistanceToNearestAtmosphereBoundary(r, mu, ray_r_mu_intersects_ground);
    float3 rayleigh_sum = 0.0;
    float3 mie_sum = 0.0;
    
    for (int i = 0; i <= SAMPLE_COUNT; ++i)
    {
        float d_i = i * dx;
        float3 rayleigh_i;
        float3 mie_i;
        ComputeSingleScatteringIntegrand(r, mu, mu_s, nu, d_i, ray_r_mu_intersects_ground, rayleigh_i, mie_i);
        float weight_i = (i == 0 || i == SAMPLE_COUNT) ? 0.5 : 1.0;
        rayleigh_sum += rayleigh_i * weight_i;
        mie_sum += mie_i * weight_i;
    }
    
    Light sun = GetMainLight();
    rayleigh = rayleigh_sum * dx * sun.color * kRayleighScattering;
    mie = mie_sum * dx * sun.color * kMieScattering;
}

void ComputeSingleScatteringTexture(float3 frag_coord, out float3 rayleigh, out float3 mie)
{
    float r, mu, mu_s, nu;
    bool ray_r_mu_intersects_ground;
    GetRMuMuSNuFromScatteringTextureFragCoord(frag_coord, r, mu, mu_s, nu, ray_r_mu_intersects_ground);
    ComputeSingleScattering(r, mu, mu_s, nu, ray_r_mu_intersects_ground, rayleigh, mie);
}

float3 GetScattering(float r, float mu, float mu_s, float nu, bool ray_r_mu_intersects_ground)
{
    float4 uvwz = GetScatteringTextureUvwzFromRMuMuSNu(r, mu, mu_s, nu, ray_r_mu_intersects_ground);
    float tex_coord_x = uvwz.x * (kScattering_NU_TextureSize - 1);
    float tex_x = floor(tex_coord_x);
    float lerp = tex_coord_x - tex_x;
    float3 uvw0 = float3((tex_x + uvwz.y) / kScattering_NU_TextureSize, uvwz.z, uvwz.w);
    float3 uvw1 = float3((tex_x + 1.0 + uvwz.y) / kScattering_NU_TextureSize, uvwz.z, uvwz.w);
    float4 scattering = SAMPLE_TEXTURE3D(_ScatteringLUT, sampler_ScatteringLUT, uvw0) * (1.0 - lerp) + SAMPLE_TEXTURE3D(_ScatteringLUT, sampler_ScatteringLUT, uvw1) * lerp;
    return scattering.rgb;
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
    
    float3 rayleigh_scattering;
    float3 mie_scattering;
    ComputeSingleScattering(r, mu, mu_s, nu, ray_r_mu_intersects_ground, rayleigh_scattering, mie_scattering);
    
    return rayleigh_scattering * RayleighPhaseFunction(nu) + mie_scattering * MiePhaseFunction(nu, _MieG);
}