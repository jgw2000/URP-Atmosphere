
Number ClampCosine(Number mu)
{
    return clamp(mu, Number(-1.0), Number(1.0));
}

Length ClampDistance(Length d)
{
    return max(d, 0.0 * m);
}

Length ClampRadius(Length r)
{
    return clamp(r, bottom_radius, top_radius);
}

Number mod(Number x, Number y)
{
    return x - y * floor(x / y);
}

float3 RadianceToLuminance(float3 rad)
{
    return mul(luminanceFromRadiance, float4(rad, 1)).rgb;
}