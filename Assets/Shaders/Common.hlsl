#pragma once

// Reference: Precomputed Atmospheric Scattering
// www.klayge.org/material/4_0/Atmospheric/Precomputed Atmospheric Scattering.pdf
static const float kUnit = 1000.0;
static const float kBottomRadius = 6360000.0 / kUnit;
static const float kTopRadius = 6420000.0 / kUnit;
static const float kRayleighScaleHeight = 8000.0 / kUnit;
static const float kMieScaleHeight = 1200.0 / kUnit;
static const float3 kRayleighScattering = float3(5.8e-6, 13.5e-6, 33.1e-6) * kUnit;
static const float3 kMieScattering = float3(21e-6, 21e-6, 21e-6) * kUnit;

static const float kTransmittance_TextureSize = 256.f;
static const float kScattering_R_TextureSize = 32.f;
static const float kScattering_MU_TextureSize = 128.f;
static const float kScattering_MU_S_TextureSize = 32.f;
static const float kScattering_NU_TextureSize = 8.f;

static const float kSunAngularRadius = 0.00935 / 2.0;
static const float3 kEarthCenter = float3(0.0, -kBottomRadius, 0.0);
static const float kmu_s_min = -0.2f;

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

float GetTextureCoordFromUnitRange(float x, float texture_size)
{
    return 0.5 / texture_size + x * (1.0 - 1.0 / texture_size);
}

float GetUnitRangeFromTextureCoord(float u, float texture_size)
{
    return (u - 0.5 / texture_size) / (1.0 - 1.0 / texture_size);
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
    return float2(GetTextureCoordFromUnitRange(x_mu, kTransmittance_TextureSize), GetTextureCoordFromUnitRange(x_r, kTransmittance_TextureSize));
}

void GetRMuFromTransmittanceTextureUv(float2 uv, out float r, out float mu)
{
    float x_mu = GetUnitRangeFromTextureCoord(uv.x, kTransmittance_TextureSize);
    float x_r = GetUnitRangeFromTextureCoord(uv.y, kTransmittance_TextureSize);
    float H = sqrt(kTopRadius * kTopRadius - kBottomRadius * kBottomRadius);
    float rho = H * x_r;
    r = sqrt(rho * rho + kBottomRadius * kBottomRadius);
    
    float d_min = kTopRadius - r;
    float d_max = rho + H;
    float d = d_min + x_mu * (d_max - d_min);
    mu = d == 0.0 ? 1.0 : (H * H - rho * rho - d * d) / (2.0 * r * d);
}

float RayleighPhaseFunction(float nu)
{
    float k = 3.0 / (16.0 * 3.1415926);
    return k * (1.0 + nu * nu);
}

float MiePhaseFunction(float nu, float g)
{
    float k = 3.0 / (8.0 * 3.1415926) * (1.0 - g * g) / (2.0 + g * g);
    return k * (1.0 + nu * nu) / pow(1.0 + g * g - 2.0 * g * nu, 1.5);
}

float4 GetScatteringTextureUvwzFromRMuMuSNu(float r, float mu, float mu_s, float nu, bool ray_r_mu_intersects_ground)
{
    float H = sqrt(kTopRadius * kTopRadius - kBottomRadius * kBottomRadius);
    float rho = sqrt(max(r * r - kBottomRadius * kBottomRadius, 0));
    float u_r = GetTextureCoordFromUnitRange(rho / H, kScattering_R_TextureSize);
    
    float r_mu = r * mu;
    float discriminant = r_mu * r_mu - r * r + kBottomRadius * kBottomRadius;
    float u_mu;
    if (ray_r_mu_intersects_ground)
    {
        float d = -r_mu - max(discriminant, 0);
        float d_min = r - kBottomRadius;
        float d_max = rho;
        u_mu = 0.5 - 0.5 * GetTextureCoordFromUnitRange(d_max == d_min ? 0.0 : (d - d_min) / (d_max - d_min), kScattering_MU_TextureSize / 2);
    }
    else
    {
        float d = -r_mu + max(discriminant + H * H, 0);
        float d_min = kTopRadius - r;
        float d_max = rho + H;
        u_mu = 0.5 + 0.5 * GetTextureCoordFromUnitRange((d - d_min) / (d_max - d_min), kScattering_MU_TextureSize / 2);
    }
    
    float d = DistanceToTopAtmosphereBoundary(kBottomRadius, mu_s);
    float d_min = kTopRadius - kBottomRadius;
    float d_max = H;
    float a = (d - d_min) / (d_max - d_min);
    float D = DistanceToTopAtmosphereBoundary(kBottomRadius, kmu_s_min);
    float A = (D - d_min) / (d_max - d_min);
    float u_mu_s = GetTextureCoordFromUnitRange(max(1.0 - a / A, 0.0) / (1.0 + a), kScattering_MU_S_TextureSize);
    
    float u_nu = (nu + 1.0) / 2.0;
    return float4(u_nu, u_mu_s, u_mu, u_r);
}

void GetRMuMuSNuFromScatteringTextureUvwz(float4 uvwz, out float r, out float mu, out float mu_s, out float nu, out bool ray_r_mu_intersects_ground)
{
    float H = sqrt(kTopRadius * kTopRadius - kBottomRadius * kBottomRadius);
    float rho = H * GetUnitRangeFromTextureCoord(uvwz.w, kScattering_R_TextureSize);
    r = sqrt(rho * rho + kBottomRadius * kBottomRadius);
    
    if (uvwz.z < 0.5)
    {
        float d_min = r - kBottomRadius;
        float d_max = rho;
        float d = d_min + (d_max - d_min) * GetUnitRangeFromTextureCoord(1.0 - 2.0 * uvwz.z, kScattering_MU_TextureSize / 2);
        mu = d == 0.0 ? -1.0 : ClampCosine(-(rho * rho + d * d) / (2.0 * r * d));
        ray_r_mu_intersects_ground = true;
    }
    else
    {
        float d_min = kTopRadius - r;
        float d_max = rho + H;
        float d = d_min + (d_max - d_min) * GetUnitRangeFromTextureCoord(2.0 * uvwz.z - 1.0, kScattering_MU_TextureSize / 2);
        mu = d == 0.0 ? 1.0 : ClampCosine((H * H - rho * rho - d * d) / (2.0 * r * d));
        ray_r_mu_intersects_ground = false;
    }

    float x_mu_s = GetUnitRangeFromTextureCoord(uvwz.y, kScattering_MU_S_TextureSize);
    float d_min = kTopRadius - kBottomRadius;
    float d_max = H;
    float D = DistanceToTopAtmosphereBoundary(kBottomRadius, kmu_s_min);
    float A = (D - d_min) / (d_max - d_min);
    float a = (A - x_mu_s * A) / (1.0 + x_mu_s * A);
    float d = d_min + min(a, A) * (d_max - d_min);
    mu_s = d == 0.0 ? 1.0 : ClampCosine((H * H - d * d) / (2.0 * kBottomRadius * d));
    
    nu = ClampCosine(uvwz.x * 2.0 - 1.0);
}

float mod(float x, float y)
{
    return x - y * floor(x / y);
}

void GetRMuMuSNuFromScatteringTextureFragCoord(float3 frag_coord, out float r, out float mu, out float mu_s, out float nu, out bool ray_r_mu_intersects_ground)
{
    const float4 SCATTERING_TEXTURE_SIZE = float4(
        kScattering_NU_TextureSize - 1,
        kScattering_MU_S_TextureSize,
        kScattering_MU_TextureSize,
        kScattering_R_TextureSize
    );

    float frag_coord_nu = floor(frag_coord.x / kScattering_MU_S_TextureSize);
    float frag_coord_mu_s = mod(frag_coord.x, kScattering_MU_S_TextureSize);
    float4 uvwz = float4(frag_coord_nu, frag_coord_mu_s, frag_coord.y, frag_coord.z) / SCATTERING_TEXTURE_SIZE;
    GetRMuMuSNuFromScatteringTextureUvwz(uvwz, r, mu, mu_s, nu, ray_r_mu_intersects_ground);
    nu = clamp(nu, mu * mu_s - sqrt((1.0 - mu * mu) * (1.0 - mu_s * mu_s)), mu * mu_s + sqrt((1.0 - mu * mu) * (1.0 - mu_s * mu_s)));
}