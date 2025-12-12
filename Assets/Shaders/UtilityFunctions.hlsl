
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

Length SafeSqrt(Area a)
{
    return sqrt(max(a, 0.0 * m2));
}

Number mod(Number x, Number y)
{
    return x - y * floor(x / y);
}