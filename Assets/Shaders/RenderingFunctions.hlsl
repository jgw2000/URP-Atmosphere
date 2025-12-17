
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
    ComputeSingleScattering(r, mu, mu_s, nu, ray_r_mu_intersects_ground, rayleigh_scattering, mie_scattering);
    
    return rayleigh_scattering * RayleighPhaseFunction(nu) + mie_scattering * MiePhaseFunction(mie_phase_function_g, nu);
}