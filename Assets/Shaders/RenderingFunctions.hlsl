
#ifdef COMBINED_SCATTERING_TEXTURES
float3 GetExtrapolatedSingleMieScattering(float4 scattering)
{
    if (scattering.r == 0.0)
        return float3(0, 0, 0);
    
    return scattering.rgb * scattering.a / scattering.r *
        (rayleigh_scattering.r / mie_scattering.r) *
        (mie_scattering / rayleigh_scattering);
}
#endif

IrradianceSpectrum GetCombinedScattering(
    Length r, Number mu, Number mu_s, Number nu,
    bool ray_r_mu_intersects_ground,
    out IrradianceSpectrum single_mie_scattering
)
{
    float4 uvwz = GetScatteringTextureUvwzFromRMuMuSNu(r, mu, mu_s, nu, ray_r_mu_intersects_ground);
    
    Number tex_coord_x = uvwz.x * Number(SCATTERING_TEXTURE_NU_SIZE - 1.0);
    Number tex_x = floor(tex_coord_x);
    Number lerp = tex_coord_x - tex_x;
    
    float3 uvw0 = float3((tex_x + uvwz.y) / Number(SCATTERING_TEXTURE_NU_SIZE), uvwz.z, uvwz.w);
    float3 uvw1 = float3((tex_x + 1.0 + uvwz.y) / Number(SCATTERING_TEXTURE_NU_SIZE), uvwz.z, uvwz.w);
    
#ifdef COMBINED_SCATTERING_TEXTURES
    float4 combined_scattering = SAMPLE_TEXTURE3D(_scattering_texture, sampler_scattering_texture, uvw0) * (1.0 - lerp) + SAMPLE_TEXTURE3D(_scattering_texture, sampler_scattering_texture, uvw1) * lerp;
    IrradianceSpectrum scattering = combined_scattering.rgb;
    single_mie_scattering = GetExtrapolatedSingleMieScattering(combined_scattering);
#else
    IrradianceSpectrum scattering = SAMPLE_TEXTURE3D(_scattering_texture, sampler_scattering_texture, uvw0).xyz * (1.0 - lerp) + SAMPLE_TEXTURE3D(_scattering_texture, sampler_scattering_texture, uvw1).xyz * lerp;
    single_mie_scattering = SAMPLE_TEXTURE3D(_single_mie_scattering_texture, sampler_single_mie_scattering_texture, uvw0).xyz * (1.0 - lerp) + SAMPLE_TEXTURE3D(_single_mie_scattering_texture, sampler_single_mie_scattering_texture, uvw1).xyz * lerp;
#endif
    
    return scattering;
}

RadianceSpectrum GetSkyRadiance(Position camera, Direction view_ray, Direction sun_direction)
{
    // Compute the distance to the top atmosphere boundary along the view ray,
    // assumiing the viewer is in space (or NaN if not intersect the atmosphere).
    Length r = length(camera);
    Length rmu = dot(camera, view_ray);
    
    // Compute the r, mu, mu_s and nu parameters needed for the texture lookups
    Number mu = rmu / r;
    Number mu_s = dot(camera, sun_direction) / r;
    Number nu = dot(view_ray, sun_direction);
    bool ray_r_mu_intersects_ground = RayIntersectsGround(r, mu);
    
    IrradianceSpectrum rayleigh_scattering;
    IrradianceSpectrum mie_scattering;
    rayleigh_scattering = GetCombinedScattering(r, mu, mu_s, nu, ray_r_mu_intersects_ground, mie_scattering);
    
    return rayleigh_scattering * RayleighPhaseFunction(nu) + mie_scattering * MiePhaseFunction(mie_phase_function_g, nu);
}