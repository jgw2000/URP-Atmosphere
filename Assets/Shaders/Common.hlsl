#pragma once

#define MAX_LOOP_ITERATIONS 30

// Reference: Precomputed Atmospheric Scattering
// www.klayge.org/material/4_0/Atmospheric/Precomputed Atmospheric Scattering.pdf
static const float Rg = 6.36e+6;
static const float Rt = 6.42e+6;
static const float Hr = 8e+3;
static const float Hm = 1.2e+3;
static const float3 betaR = float3(5.8e-6, 13.5e-6, 33.1e-6); // rayleigh scattering coefficients
static const float3 betaM = float3(21e-6, 21e-6, 21e-6); // mie scattering coefficients

float RaySphere(float3 sphereCenter, float sphereRadius, float3 rayOrigin, float3 rayDir)
{
    float3 offset = rayOrigin - sphereCenter;
    float a = 1; // Set to dot(rayDir, rayDir) if rayDir might not be normalized
    float b = 2 * dot(offset, rayDir);
    float c = dot(offset, offset) - sphereRadius * sphereRadius;
    float d = b * b - 4 * a * c;
    float s = sqrt(d);
    return (-b + s) / (2 * a); // here we ensure that rayOrigin is inside sphere
}

float3 DensityAtPoint(float3 sphereCenter, float3 position)
{
    float h = length(position - sphereCenter) - Rg;
    float rayleighDensity = exp(-h / Hr);
    float mieDensity = exp(-h / Hm);
    
    return float3(rayleighDensity, mieDensity, 0);
}