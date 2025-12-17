
void ComputeSingleScatteringIntegrand(
    Length r, Number mu, Number mu_s, Number nu, Length d, bool ray_r_mu_intersects_ground,
    out DimensionlessSpectrum rayleigh, out DimensionlessSpectrum mie    
)
{
    Length r_d = ClampRadius(sqrt(d * d + 2.0 * r * mu * d + r * r));
    Number mu_s_d = ClampCosine((r * mu_s + d * nu) / r_d);
    
    DimensionlessSpectrum transmittance =
        GetTransmittance(r, mu, d, ray_r_mu_intersects_ground) *
        GetTransmittanceToSun(r_d, mu_s_d);
    
    rayleigh = transmittance * GetProfileDensity(RayleighDensity(), r_d - bottom_radius);
    mie = transmittance * GetProfileDensity(MieDensity(), r_d - bottom_radius);
}

Length DistanceToNearestAtmosphereBoundary(Length r, Number mu, bool ray_r_mu_intersects_ground)
{
    if (ray_r_mu_intersects_ground)
        return DistanceToBottomAtmosphereBoundary(r, mu);
    else
        return DistanceToTopAtmosphereBoundary(r, mu);
}

void ComputeSingleScattering(
    Length r, Number mu, Number mu_s, Number nu, bool ray_r_mu_intersects_ground,
    out IrradianceSpectrum rayleigh, out IrradianceSpectrum mie
)
{
    const int SAMPLE_COUNT = 50;
    Length dx = DistanceToNearestAtmosphereBoundary(r, mu, ray_r_mu_intersects_ground) / Number(SAMPLE_COUNT);
    
    DimensionlessSpectrum rayleigh_sum = DimensionlessSpectrum(0, 0, 0);
    DimensionlessSpectrum mie_sum = DimensionlessSpectrum(0, 0, 0);
    
    for (int i = 0; i <= SAMPLE_COUNT; ++i)
    {
        Length d_i = Number(i) * dx;
        DimensionlessSpectrum rayleigh_i;
        DimensionlessSpectrum mie_i;
        ComputeSingleScatteringIntegrand(r, mu, mu_s, nu, d_i, ray_r_mu_intersects_ground, rayleigh_i, mie_i);
        
        // Sample weight
        Number weight_i = (i == 0 || i == SAMPLE_COUNT) ? 0.5 : 1.0;
        rayleigh_sum += rayleigh_i * weight_i;
        mie_sum += mie_i * weight_i;
    }

    rayleigh = rayleigh_sum * dx * solar_irradiance * rayleigh_scattering;
    mie = mie_sum * dx * solar_irradiance * mie_scattering;
}

InverseSolidAngle RayleighPhaseFunction(Number nu)
{
    InverseSolidAngle k = 3.0 / (16.0 * PI * sr);
    return k * (1.0 + nu * nu);
}

InverseSolidAngle MiePhaseFunction(Number g, Number nu)
{
    InverseSolidAngle k = 3.0 / (8.0 * PI * sr) * (1.0 - g * g) / (2.0 + g * g);
    return k * (1.0 + nu * nu) / pow(abs(1.0 + g * g - 2.0 * g * nu), 1.5);
}