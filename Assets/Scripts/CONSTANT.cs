using UnityEngine;

public enum LUMINANCE
{
    // Render the spectral radiance at kLambdaR, kLambdaG, kLambdaB.
    NONE,
    // Render the sRGB luminance, using an approximate (on the fly) conversion
    // from 3 spectral radiance values only (see section 14.3 in <a href=
    // "https://arxiv.org/pdf/1612.04336.pdf">A Qualitative and Quantitative
    //  Evaluation of 8 Clear Sky Models</a>).
    APPROXIMATE,
    // Render the sRGB luminance, precomputed from 15 spectral radiance values
    // (see section 4.4 in <a href=
    // "http://www.oskee.wz.cz/stranka/uploads/SCCG10ElekKmoch.pdf">Real-time
    //  Spectral Scattering in Large-scale Natural Participating Media</a>).
    PRECOMPUTED
}

public static class CONSTANTS
{
    public static readonly int NUM_THREADS = 8;

    public static readonly int TRANSMITTANCE_WIDTH = 256;
    public static readonly int TRANSMITTANCE_HEIGHT = 64;
    public static readonly int TRANSMITTANCE_CHANNELS = 3;
    public static readonly int TRANSMITTANCE_SIZE = TRANSMITTANCE_WIDTH * TRANSMITTANCE_HEIGHT;

    public static readonly int SCATTERING_R = 32;
    public static readonly int SCATTERING_MU = 128;
    public static readonly int SCATTERING_MU_S = 32;
    public static readonly int SCATTERING_NU = 8;

    public static readonly int SCATTERING_WIDTH = SCATTERING_NU * SCATTERING_MU_S;
    public static readonly int SCATTERING_HEIGHT = SCATTERING_MU;
    public static readonly int SCATTERING_DEPTH = SCATTERING_R;
    public static readonly int SCATTERING_CHANNELS = 4;
    public static readonly int SCATTERING_SIZE = SCATTERING_WIDTH * SCATTERING_HEIGHT * SCATTERING_DEPTH;
}
