
IrradianceSpectrum ComputeDirectIrradiance(Length r, Number mu_s)
{
    Number alpha_s = sun_angular_radius / rad;
    
    // Approximate average of the cosine factor mu_s over the visible fraction of
    // the sun disc
    Number average_cosine_factor =
        mu_s < -alpha_s ? 0.0 : (mu_s > alpha_s ? mu_s :
            (mu_s + alpha_s) * (mu_s + alpha_s) / (4.0 * alpha_s));
    
    return solar_irradiance * GetTransmittanceToTopAtmosphereBoundary(r, mu_s) * average_cosine_factor;
}

float2 GetIrradianceTextureUvFromRMuS(Length r, Number mu_s)
{
    Number x_r = (r - bottom_radius) / (top_radius - bottom_radius);
    Number x_mu_s = mu_s * 0.5 + 0.5;
    return float2(GetTextureCoordFromUnitRange(x_mu_s, IRRADIANCE_TEXTURE_WIDTH),
                  GetTextureCoordFromUnitRange(x_r, IRRADIANCE_TEXTURE_HEIGHT));
}

void GetRMuSFromIrradianceTextureUv(float2 uv, out Length r, out Number mu_s)
{
    Number x_mu_s = GetUnitRangeFromTextureCoord(uv.x, IRRADIANCE_TEXTURE_WIDTH);
    Number x_r = GetUnitRangeFromTextureCoord(uv.y, IRRADIANCE_TEXTURE_HEIGHT);
    r = bottom_radius + x_r * (top_radius - bottom_radius);
    mu_s = ClampCosine(2.0 * x_mu_s - 1.0);
}

IrradianceSpectrum ComputeDirectIrradianceTexture(float2 gl_frag_coord)
{
    Length r;
    Number mu_s;
    
    const float2 IRRADIANCE_TEXTURE_SIZE = float2(IRRADIANCE_TEXTURE_WIDTH, IRRADIANCE_TEXTURE_HEIGHT);
    
    GetRMuSFromIrradianceTextureUv(gl_frag_coord / IRRADIANCE_TEXTURE_SIZE, r, mu_s);
    return ComputeDirectIrradiance(r, mu_s);
}