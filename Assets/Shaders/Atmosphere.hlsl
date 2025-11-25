#pragma once

#include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/Lighting.hlsl"
#include "Common.hlsl"

static const int _InScatteringPoints = 20;

CBUFFER_START(UnityPerMaterial)
float _MieG;
CBUFFER_END

TEXTURE2D(_TransmittanceLUT);
SAMPLER(sampler_TransmittanceLUT);

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