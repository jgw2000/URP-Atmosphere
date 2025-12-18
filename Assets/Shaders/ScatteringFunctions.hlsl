
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

float4 GetScatteringTextureUvwzFromRMuMuSNu(
    Length r, Number mu, Number mu_s, Number nu,
    bool ray_r_mu_intersects_ground
)
{
    // Distance to top atmosphere boundary for a horizontal ray at ground level
    Length H = sqrt(top_radius * top_radius - bottom_radius * bottom_radius);
    
    // Distance to the horizon
    Length rho = SafeSqrt(r * r - bottom_radius * bottom_radius);
    Number u_r = GetTextureCoordFromUnitRange(rho / H, SCATTERING_TEXTURE_R_SIZE);
    
    // Discriminant of the quadratic equation for the intersections of the ray
    // (r,mu) with the ground
    Length r_mu = r * mu;
    Area discriminant = r_mu * r_mu - r * r + bottom_radius * bottom_radius;
    Number u_mu;
    
    if (ray_r_mu_intersects_ground)
    {
        // Distance to the ground for the ray (r,mu), and its minimum and maximum
        // values over all mu - obtained for (r,-l) and (r,mu_horizon)
        Length d = -r_mu - SafeSqrt(discriminant);
        Length d_min = r - bottom_radius;
        Length d_max = rho;
        u_mu = 0.5 - 0.5 * GetTextureCoordFromUnitRange(d_max == d_min ? 0.0 :
            (d - d_min) / (d_max - d_min), SCATTERING_TEXTURE_MU_SIZE / 2.0);
    }
    else
    {
        // Distance to the top atmosphere boundary for the ray (r,mu), and its
        // minimum and maximum values over all mu - obtained for (r,l) and
        // (r,mu_horizon)
        Length d = -r_mu + SafeSqrt(discriminant + H * H);
        Length d_min = top_radius - r;
        Length d_max = rho + H;
        u_mu = 0.5 + 0.5 * GetTextureCoordFromUnitRange(
            (d - d_min) / (d_max - d_min), SCATTERING_TEXTURE_MU_SIZE / 2.0);
    }
    
    Length d = DistanceToTopAtmosphereBoundary(bottom_radius, mu_s);
    Length d_min = top_radius - bottom_radius;
    Length d_max = H;
    Number a = (d - d_min) / (d_max - d_min);
    // Length D = DistanceToTopAtmosphereBoundary(bottom_radius, mu_s_min);
    // Number A = (D - d_min) / (d_max - d_min);
    Number A = -2.0 * mu_s_min * bottom_radius / (d_max - d_min);
    
    // An ad-hoc function equal to 0 for mu_s = mu_s_min (because then d = D and
    // thus a = A), equal to 1 for mu_s = 1 (because then d = d_min and thus
    // a = 0), and with a large slope around mu_s = 0, to get more texture
    // samples near the horizon
    Number u_mu_s = GetTextureCoordFromUnitRange(max(1.0 - a / A, 0.0) / (1.0 + a), SCATTERING_TEXTURE_MU_S_SIZE);
    
    Number u_nu = (nu + 1.0) / 2.0;
    return float4(u_nu, u_mu_s, u_mu, u_r);
}

void GetRMuMuSNuFromScatteringTextureUvwz(
    float4 uvwz, out Length r, out Number mu, out Number mu_s,
    out Number nu, out bool ray_r_mu_intersects_ground
)
{
    // Distance to top atmosphere boundary for a horizontal ray at ground level
    Length H = sqrt(top_radius * top_radius - bottom_radius * bottom_radius);
    // Distance to the horizon
    Length rho = H * GetUnitRangeFromTextureCoord(uvwz.w, SCATTERING_TEXTURE_R_SIZE);
    r = sqrt(rho * rho + bottom_radius * bottom_radius);
    
    if (uvwz.z < 0.5)
    {
        // Distance to the ground for the ray (r,mu), and its minimum and maximum
        // values over all mu - obtained for (r,-l) and (r,mu_horizon)
        Length d_min = r - bottom_radius;
        Length d_max = rho;
        Length d = d_min + (d_max - d_min) * GetUnitRangeFromTextureCoord(1.0 - 2.0 * uvwz.z, SCATTERING_TEXTURE_MU_SIZE / 2.0);
        mu = d == 0.0 * m ? Number(-1.0) : ClampCosine(-(rho * rho + d * d) / (2.0 * r * d));
        ray_r_mu_intersects_ground = true;
    }
    else
    {
        // Distance to the top atmosphere boundary for the ray (r,mu), and its
        // minimum and maximum values over all mu - obtained for (r,l) and
        // (r,mu_horizon)
        Length d_min = top_radius - r;
        Length d_max = rho + H;
        Length d = d_min + (d_max - d_min) * GetUnitRangeFromTextureCoord(2.0 * uvwz.z - 1.0, SCATTERING_TEXTURE_MU_SIZE / 2.0);
        mu = d == 0.0 * m ? Number(1.0) : ClampCosine((H * H - rho * rho - d * d) / (2.0 * r * d));
        ray_r_mu_intersects_ground = false;
    }
    
    Number x_mu_s = GetUnitRangeFromTextureCoord(uvwz.y, SCATTERING_TEXTURE_MU_S_SIZE);
    Length d_min = top_radius - bottom_radius;
    Length d_max = H;
    // Length D = DistanceToTopAtmosphereBoundary(bottom_radius, mu_s_min);
    // Number A = (D - d_min) / (d_max - d_min);
    Number A = -2.0 * mu_s_min * bottom_radius / (d_max - d_min);
    Number a = (A - x_mu_s * A) / (1.0 + x_mu_s * A);
    Length d = d_min + min(a, A) * (d_max - d_min);
    mu_s = d == 0.0 * m ? Number(1.0) : ClampCosine((H * H - d * d) / (2.0 * bottom_radius * d));
    
    nu = ClampCosine(uvwz.x * 2.0 - 1.0);
}

/*
 * We assumed above that we have 4D textures, which is not the case in practice.
 * We therefore need a further mapping, between 3D and 4D texture coordinates. The
 * function below expands a 3D texel coordinate into a 4D texture coordinate, and
 * then to (r,mu,mu_s,nu) parameters. It does so by "unpacking" two texel
 * coordinates from the x texel coordinate. Note also how we clamp the nu
 * parameter at the end. This is because nu is not a fully independent variable:
 * its range of values depends on mu and mu_s (this can be seen by computing
 * mu, mu_s and nu from the cartesian coordinates of the zenith, view and
 * sun unit direction vectors), and the previous functions implicitely assume this
 * (their assertions can break if this constraint is not respected).
 */

void GetRMuMuSNuFromScatteringTextureFragCoord(
    float3 gl_frag_coord,
    out Length r, out Number mu, out Number mu_s, out Number nu,
    out bool ray_r_mu_intersects_ground
)
{
    const float4 SCATTERING_TEXTURE_SIZE = float4(
        SCATTERING_TEXTURE_NU_SIZE - 1,
        SCATTERING_TEXTURE_MU_S_SIZE,
        SCATTERING_TEXTURE_MU_SIZE,
        SCATTERING_TEXTURE_R_SIZE
    );
    
    Number frag_coord_nu = floor(gl_frag_coord.x / Number(SCATTERING_TEXTURE_MU_S_SIZE));
    Number frag_coord_mu_s = mod(gl_frag_coord.x, Number(SCATTERING_TEXTURE_MU_S_SIZE));
    
    float4 uvwz = float4(frag_coord_nu, frag_coord_mu_s, gl_frag_coord.y, gl_frag_coord.z) / SCATTERING_TEXTURE_SIZE;
    
    GetRMuMuSNuFromScatteringTextureUvwz(uvwz, r, mu, mu_s, nu, ray_r_mu_intersects_ground);
    
    // Clamp nu to its valid range of values, given mu and mu_s
    nu = clamp(nu, mu * mu_s - sqrt((1.0 - mu * mu) * (1.0 - mu_s * mu_s)),
        mu * mu_s + sqrt((1.0 - mu * mu) * (1.0 - mu_s * mu_s)));
}

void ComputeSingleScatteringTexture(
    float3 gl_frag_coord,
    out IrradianceSpectrum rayleigh, out IrradianceSpectrum mie
)
{
    Length r;
    Number mu;
    Number mu_s;
    Number nu;
    bool ray_r_mu_intersects_ground;
    GetRMuMuSNuFromScatteringTextureFragCoord(gl_frag_coord, r, mu, mu_s, nu, ray_r_mu_intersects_ground);
    ComputeSingleScattering(r, mu, mu_s, nu, ray_r_mu_intersects_ground, rayleigh, mie);
}

/*
 * Lookup
 * 
 * With the help of the above precomputed texture, we can now get the scattering
 * between a point and the nearest atmosphere boundary with two texture lookups (we
 * need two 3D texture lookups to emulate a single 4D texture lookup with
 * quadrilinear interpolation; the 3D texture coordinates are computed using the
 * inverse of the 3D-4D mapping defined in
 * GetRMuMuSNuFromScatteringTextureFragCoord):
 */
AbstractSpectrum GetScattering(
    Length r, Number mu, Number mu_s, Number nu,
    bool ray_r_mu_intersects_ground
)
{
    float4 uvwz = GetScatteringTextureUvwzFromRMuMuSNu(r, mu, mu_s, nu, ray_r_mu_intersects_ground);
    
    Number tex_coord_x = uvwz.x * Number(SCATTERING_TEXTURE_NU_SIZE - 1);
    Number tex_x = floor(tex_coord_x);
    Number lerp = tex_coord_x - tex_x;
    
    float3 uvw0 = float3((tex_x + uvwz.y) / Number(SCATTERING_TEXTURE_NU_SIZE), uvwz.z, uvwz.w);
    float3 uvw1 = float3((tex_x + 1.0 + uvwz.y) / Number(SCATTERING_TEXTURE_NU_SIZE), uvwz.z, uvwz.w);
    
    return SAMPLE_TEXTURE3D(_scattering_texture, sampler_scattering_texture, uvw0).rgb * (1.0 - lerp) +
           SAMPLE_TEXTURE3D(_scattering_texture, sampler_scattering_texture, uvw1).rgb * lerp;
}