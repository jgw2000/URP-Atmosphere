#pragma once

#define MAX_LOOP_ITERATIONS 30

// Reference: Precomputed Atmospheric Scattering
// www.klayge.org/material/4_0/Atmospheric/Precomputed Atmospheric Scattering.pdf
static const float kUnit = 1000.0;
static const float kBottomRadius = 6360000.0 / kUnit;
static const float kTopRadius = 6420000.0 / kUnit;
static const float kRayleighScaleHeight = 8000.0 / kUnit;
static const float kMieScaleHeight = 1200.0 / kUnit;
static const float3 kRayleighScattering = float3(5.8e-6, 13.5e-6, 33.1e-6) * kUnit;
static const float3 kMieScattering = float3(21e-6, 21e-6, 21e-6) * kUnit;

static const float kTexelSize = 1.0 / 256.f;
static const float kSunAngularRadius = 0.00935 / 2.0;
static const float3 kEarthCenter = float3(0.0, -kBottomRadius, 0.0);

static const float Rg = 6.36;
static const float Rt = 6.42;
static const float Hr = 8e-3;
static const float Hm = 1.2e-3;
static const float3 betaR = float3(5.8, 13.5, 33.1); // rayleigh scattering coefficients
static const float3 betaM = float3(21, 21, 21); // mie scattering coefficients

float ClampCosine(float mu)
{
    return clamp(mu, -1.0, 1.0);
}

float ClampDistance(float d)
{
    return max(d, 0.0);
}

float ClampRadius(float r)
{
    return clamp(r, kBottomRadius, kTopRadius);
}

float DistanceToTopAtmosphereBoundary(float r, float mu)
{
    float discriminant = r * r * (mu * mu - 1.0) + kTopRadius * kTopRadius;
    return ClampDistance(-r * mu + sqrt(max(discriminant, 0)));
}

float DistanceToBottomAtmosphereBoundary(float r, float mu)
{
    float discriminant = r * r * (mu * mu - 1.0) + kBottomRadius * kBottomRadius;
    return ClampDistance(-r * mu - sqrt(max(discriminant, 0)));
}

bool RayIntersectsGround(float r, float mu)
{
    float discriminant = r * r * (mu * mu - 1.0) + kBottomRadius * kBottomRadius;
    return mu < 0.0 && discriminant >= 0.0;
}

float GetLayerDensity(float scaledHeight, float altitude)
{
    float density = exp(-altitude / scaledHeight);
    return clamp(density, 0.0, 1.0);
}

float ComputeOpticalLengthToTopAtmosphereBoundary(float r, float mu)
{
    const int SAMPLE_COUNT = 500;
    float dx = DistanceToTopAtmosphereBoundary(r, mu) / SAMPLE_COUNT;
    float result = 0.0;
    
    for (int i = 0; i <= SAMPLE_COUNT; ++i)
    {
        float d_i = i * dx;
        float r_i = sqrt(d_i * d_i + 2.0 * r * mu * d_i + r * r);
        float y_i = GetLayerDensity(kRayleighScaleHeight, r_i - kBottomRadius);
        float weight_i = i == 0 || i == SAMPLE_COUNT ? 0.5 : 1.0;
        result += y_i * weight_i * dx;
    }
    
    return result;
}

float3 ComputeTransmittanceToTopAtmosphereBoundary(float r, float mu)
{
    return exp(-kRayleighScattering * ComputeOpticalLengthToTopAtmosphereBoundary(r, mu));
}

float GetTextureCoordFromUnitRange(float x)
{
    return 0.5 * kTexelSize + x * (1.0 - 1.0 * kTexelSize);
}

float GetUnitRangeFromTextureCoord(float u)
{
    return (u - 0.5 * kTexelSize) / (1.0 - 1.0 * kTexelSize);
}

float2 GetTransmittanceTextureUvFromRMu(float r, float mu)
{
    float H = sqrt(kTopRadius * kTopRadius - kBottomRadius * kBottomRadius);
    float rho = sqrt(max(r * r - kBottomRadius * kBottomRadius, 0));
    float d = DistanceToTopAtmosphereBoundary(r, mu);
    float d_min = kTopRadius - r;
    float d_max = rho + H;
    float x_mu = (d - d_min) / (d_max - d_min);
    float x_r = rho / H;
    return float2(GetTextureCoordFromUnitRange(x_mu), GetTextureCoordFromUnitRange(x_r));
}

float2 GetRMuFromTransmittanceTextureUv(float2 uv)
{
    float x_mu = GetUnitRangeFromTextureCoord(uv.x);
    float x_r = GetUnitRangeFromTextureCoord(uv.y);
    float H = sqrt(kTopRadius * kTopRadius - kBottomRadius * kBottomRadius);
    float rho = H * x_r;
    float r = sqrt(rho * rho + kBottomRadius * kBottomRadius);
    float d_min = kTopRadius - r;
    float d_max = rho + H;
    float d = d_min + x_mu * (d_max - d_min);
    float mu = d == 0.0 ? 1.0 : (H * H - rho * rho - d * d) / (2.0 * r * d);
    return float2(r, ClampCosine(mu));
}

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