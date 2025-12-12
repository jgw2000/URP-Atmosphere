
// Mapping x in [0, 1] to texture coordinates u in [0.5/n, 1-0.5/n]
Number GetTextureCoordFromUnitRange(Number x, int texture_size)
{
    return 0.5 / Number(texture_size) + x * (1.0 - 1.0 / Number(texture_size));
}

// Mapping texture coordinates u in [0.5/n, 1-0.5/n] to x in [0, 1]
Number GetUnitRangeFromTextureCoord(Number u, int texture_size)
{
    return (u - 0.5 / Number(texture_size)) / (1.0 - 1.0 / Number(texture_size));
}

void GetRMuFromTransmittanceTextureUv(float2 uv, out Length r, out Number mu)
{
    Number x_mu = GetUnitRangeFromTextureCoord(uv.x, TRANSMITTANCE_TEXTURE_WIDTH);
    Number x_r = GetUnitRangeFromTextureCoord(uv.y, TRANSMITTANCE_TEXTURE_HEIGHT);
    // Distance to top atmosphere boundary for a horizontal ray at ground level.
    Length H = sqrt(top_radius * top_radius - bottom_radius * bottom_radius);
    // Distance to the horizon, from which we can compute r:
    Length rho = H * x_r;
    r = sqrt(rho * rho + bottom_radius * bottom_radius);
    // Distance to the top atmosphere boundary for the ray (r,mu), and its minimum
    // and maximum values over all mu - obtained for (r,l) and (r,mu_horizon) -
    // from which we can recover mu:
    Length d_min = top_radius - r;
    Length d_max = rho + H;
    Length d = d_min + x_mu * (d_max - d_min);
    mu = d == 0.0 * m ? Number(1.0) : (H * H - rho * rho - d * d) / (2.0 * r * d);
    mu = ClampCosine(mu);
}

Length DistanceToTopAtmosphereBoundary(Length r, Number mu)
{
    Area discriminant = r * r * (mu * mu - 1.0) + top_radius * top_radius;
    return ClampDistance(-r * mu + SafeSqrt(discriminant));
}

Length DistanceToBottomAtmosphereBoundary(Length r, Number mu)
{
    Area discriminant = r * r * (mu * mu - 1.0) + bottom_radius * bottom_radius;
    return ClampDistance(-r * mu - SafeSqrt(discriminant));
}

Number GetLayerDensity(DensityProfileLayer layer, Length altitude)
{
    Number density = layer.exp_term * exp(layer.exp_scale * altitude) + layer.linear_term * altitude + layer.constant_term;
    return clamp(density, Number(0.0), Number(1.0));
}

Number GetProfileDensity(DensityProfile profile, Length altitude)
{
    return altitude < profile.layers[0].width ?
        GetLayerDensity(profile.layers[0], altitude) :
        GetLayerDensity(profile.layers[1], altitude);
}

Length ComputeOpticalLengthToTopAtmosphereBoundary(DensityProfile profile, Length r, Number mu)
{
    const int SAMPLE_COUNT = 500;
    Length dx = DistanceToTopAtmosphereBoundary(r, mu) / Number(SAMPLE_COUNT);
    Length result = 0.0 * m;
    
    for (int i = 0; i <= SAMPLE_COUNT; ++i)
    {
        Length d_i = Number(i) * dx;
        // Distance between the current sample point and the planet center
        Length r_i = sqrt(d_i * d_i + 2.0 * r * mu * d_i + r * r);
        // Number density at the current sample point
        Number y_i = GetProfileDensity(profile, r_i - bottom_radius);
        Number weight_i = i == 0 || i == SAMPLE_COUNT ? 0.5 : 1.0;
        result += y_i * weight_i * dx;
    }
    
    return result;
}

// Assume the segment does not intersect the ground
DimensionlessSpectrum ComputeTransmittanceToTopAtmosphereBoundary(Length r, Number mu)
{
    DimensionlessSpectrum density = 0;
    density += rayleigh_scattering * ComputeOpticalLengthToTopAtmosphereBoundary(RayleighDensity(), r, mu);
    density += mie_scattering * ComputeOpticalLengthToTopAtmosphereBoundary(MieDensity(), r, mu);
    density += absorption_extinction * ComputeOpticalLengthToTopAtmosphereBoundary(AbsorptionDensity(), r, mu);
    
    return exp(-density);
}

// Precompute a texel of the transmittance texture
DimensionlessSpectrum ComputeTransmittanceToTopAtmosphereBoundaryTexture(float2 gl_frag_coord)
{
    const float2 TRANSMITTANCE_TEXTURE_SIZE = float2(TRANSMITTANCE_TEXTURE_WIDTH, TRANSMITTANCE_TEXTURE_HEIGHT);
    Length r;
    Number mu;
    GetRMuFromTransmittanceTextureUv(gl_frag_coord / TRANSMITTANCE_TEXTURE_SIZE, r, mu);
    return ComputeTransmittanceToTopAtmosphereBoundary(r, mu);
}