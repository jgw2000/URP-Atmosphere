/**
 * Copyright (c) 2017 Eric Bruneton
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holders nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 */

/* This GLSL file defines the physical types and constants which are used in the
 * main functions of our atmosphere model.
 *
 * Physical quantities
 *
 * The physical quantities we need for our atmosphere model are radiometric and photometric quantities.
 * 
 * https://en.wikipedia.org/wiki/Radiometry
 * https://en.wikipedia.org/wiki/Photometry_(optics)
 *
 * We start with six base quantities: length, wavelength, angle, solid angle,
 * power and luminous power (wavelength is also a length, but we distinguish the
 * two for increased clarity).
 *
 */

#define Length float
#define Wavelength float
#define Angle float
#define SolidAngle float
#define Power float
#define LuminousPower float

/*
 * From this we "derive" the irradiance, radiance, spectral irradiance,
 * spectral radiance, luminance, etc, as well pure numbers, area, volume.
 */

#define Number float
#define InverseLength float
#define Area float
#define Volume float
#define NumberDensity float
#define Irradiance float
#define Radiance float
#define SpectralPower float
#define SpectralIrradiance float
#define SpectralRadiance float
#define SpectralRadianceDensity float
#define ScatteringCoefficient float
#define InverseSolidAngle float
#define LuminousIntensity float
#define Luminance float
#define Illuminance float

/*
 * We  also need vectors of physical quantities, mostly to represent functions
 * depending on the wavelength. In this case the vector elements correspond to
 * values of a function at some predefined wavelengths.
 */

// A generic function from Wavelength to some other type.
#define AbstractSpectrum float3
// A function from Wavelength to Number.
#define DimensionlessSpectrum float3
// A function from Wavelength to SpectralPower.
#define PowerSpectrum float3
// A function from Wavelength to SpectralIrradiance.
#define IrradianceSpectrum float3
// A function from Wavelength to SpectralRadiance.
#define RadianceSpectrum float3
// A function from Wavelength to SpectralRadianceDensity.
#define RadianceDensitySpectrum float3
// A function from Wavelength to ScaterringCoefficient.
#define ScatteringSpectrum float3

// A position in 3D (3 length values).
#define Position float3
// A unit direction vector in 3D (3 unitless values).
#define Direction float3
// A vector of 3 luminance values.
#define Luminance3 float3
// A vector of 3 illuminance values.
#define Illuminance3 float3

/*
 * Finally, we also need precomputed textures containing physical quantities in each texel.
 */

#ifdef COMPUTE_SHADER

#define TransmittanceTexture Texture2D<float4>
#define AbstractScatteringTexture Texture3D<float4>
#define ReducedScatteringTexture Texture3D<float4>
#define ScatteringTexture Texture3D<float4>
#define ScatteringDensityTexture Texture3D<float4>
#define IrradianceTexture Texture2D<float4>

#else

#define TransmittanceTexture sampler2D
#define AbstractScatteringTexture sampler3D
#define ReducedScatteringTexture sampler3D
#define ScatteringTexture sampler3D
#define ScatteringDensityTexture sampler3D
#define IrradianceTexture sampler2D

#endif

/*
 * Physical units
 *
 * We can then define the units for our six base physical quantities:
 * meter (m), nanometer (nm), radian (rad), steradian (sr), watt (watt) and lumen (lm):
 */

static const Length m = 1.0;
static const Wavelength nm = 1.0;
static const Angle rad = 1.0;
static const SolidAngle sr = 1.0;
static const Power watt = 1.0;
static const LuminousPower lm = 1.0;

/*
 * From which we can derive the units for some derived physical quantities,
 * as well as some derived units (kilometer km, kilocandela kcd, degree deg):
 */

static const float PI = 3.14159265358979323846;

static const Length km = 1000.0 * m;
static const Area m2 = m * m;
static const Volume m3 = m * m * m;
static const Angle pi = PI * rad;
static const Angle deg = pi / 180.0;
static const Irradiance watt_per_square_meter = watt / m2;
static const Radiance watt_per_square_meter_per_sr = watt / (m2 * sr);
static const SpectralIrradiance watt_per_square_meter_per_nm = watt / (m2 * nm);
static const SpectralRadiance watt_per_square_meter_per_sr_per_nm = watt / (m2 * sr * nm);
static const SpectralRadianceDensity watt_per_cubic_meter_per_sr_per_nm = watt / (m3 * sr * nm);
static const LuminousIntensity cd = lm / sr;
static const LuminousIntensity kcd = 1000.0 * cd;
static const Luminance cd_per_square_meter = cd / m2;
static const Luminance kcd_per_square_meter = kcd / m2;

/*
 * Atmosphere parameters
 *
 * Using the above types, we can now define the parameters of our atmosphere
 * model. We start with the definition of density profiles, which are needed for
 * parameters that depend on the altitude:
 */

int TRANSMITTANCE_TEXTURE_WIDTH;
int TRANSMITTANCE_TEXTURE_HEIGHT;

// Rayleigh
ScatteringSpectrum rayleigh_scattering;

Length rayleigh_width;
Number rayleigh_exp_term;
InverseLength rayleigh_exp_scale;
InverseLength rayleigh_linear_term;
Number rayleigh_constant_term;

// Mie
ScatteringSpectrum mie_scattering;
ScatteringSpectrum mie_extinction;

Length mie_width;
Number mie_exp_term;
InverseLength mie_exp_scale;
InverseLength mie_linear_term;
Number mie_constant_term;

// Ozone
ScatteringSpectrum absorption_extinction;

Length absorption0_width;
Number absorption0_exp_term;
InverseLength absorption0_exp_scale;
InverseLength absorption0_linear_term;
Number absorption0_constant_term;

Length absorption1_width;
Number absorption1_exp_term;
InverseLength absorption1_exp_scale;
InverseLength absorption1_linear_term;
Number absorption1_constant_term;

// Sun

Length bottom_radius;
Length top_radius;

struct DensityProfileLayer
{
    Length width;
    Number exp_term;
    InverseLength exp_scale;
    InverseLength linear_term;
    Number constant_term;
};

struct DensityProfile
{
    DensityProfileLayer layers[2];
};

DensityProfile RayleighDensity()
{
    DensityProfileLayer layer;
    
    layer.width = rayleigh_width;
    layer.exp_term = rayleigh_exp_term;
    layer.exp_scale = rayleigh_exp_scale;
    layer.linear_term = rayleigh_linear_term;
    layer.constant_term = rayleigh_constant_term;
    
    DensityProfile profile;
    profile.layers[0] = layer;
    profile.layers[1] = layer;
    
    return profile;
}

DensityProfile MieDensity()
{
    DensityProfileLayer layer;
    
    layer.width = mie_width;
    layer.exp_term = mie_exp_term;
    layer.exp_scale = mie_exp_scale;
    layer.linear_term = mie_linear_term;
    layer.constant_term = mie_constant_term;
    
    DensityProfile profile;
    profile.layers[0] = layer;
    profile.layers[1] = layer;
    
    return profile;
}

DensityProfile AbsorptionDensity()
{
    DensityProfileLayer layer0;
    
    layer0.width = absorption0_width;
    layer0.exp_term = absorption0_exp_term;
    layer0.exp_scale = absorption0_exp_scale;
    layer0.linear_term = absorption0_linear_term;
    layer0.constant_term = absorption0_constant_term;
    
    DensityProfileLayer layer1;
    
    layer1.width = absorption1_width;
    layer1.exp_term = absorption1_exp_term;
    layer1.exp_scale = absorption1_exp_scale;
    layer1.linear_term = absorption1_linear_term;
    layer1.constant_term = absorption1_constant_term;
    
    DensityProfile profile;
    profile.layers[0] = layer0;
    profile.layers[1] = layer1;
    
    return profile;
}