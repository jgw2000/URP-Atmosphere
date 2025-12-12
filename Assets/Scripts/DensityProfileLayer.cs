using UnityEngine;

// ¦Ñ(h) = clamp(exp_term * e^(exp_scale * h) + linear_term * h + constant_term, 0, 1)
public struct DensityProfileLayer
{
    public string Name { get; private set; }

    public double Width { get; private set; }

    public double ExpTerm { get; private set; }

    public double ExpScale { get; private set; }

    public double LinearTerm { get; private set; }

    public double ConstantTerm { get; private set; }

    public DensityProfileLayer(string name, double width, double exp_term, double exp_scale, double linear_term, double constant_term)
    {
        Name = name;
        Width = width;
        ExpTerm = exp_term;
        ExpScale = exp_scale;
        LinearTerm = linear_term;
        ConstantTerm = constant_term;
    }
}
