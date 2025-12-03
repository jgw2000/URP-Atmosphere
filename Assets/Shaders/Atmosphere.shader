Shader "Custom/Atmosphere"
{
    Properties
    {
        [NoScaleOffset] _TransmittanceLUT("Transmittance LUT", 2D) = "black" {}
        _MieG("Mie G", Range(-1, 1)) = 0.9
    }

    SubShader
    {
        Tags { "RenderType"="Opaque" "RenderPipeline"="UniversalPipeline" }
        LOD 100
        ZWrite Off Cull Off

        Pass
        {
            HLSLPROGRAM
            #include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/Core.hlsl"
            #include "Packages/com.unity.render-pipelines.core/Runtime/Utilities/Blit.hlsl"
            #include "Atmosphere.hlsl"

            #pragma vertex Vert
            #pragma fragment frag

            float4 frag(Varyings input) : SV_Target
            {
                float4 viewPos = mul(unity_CameraInvProjection, float4(input.texcoord * 2 - 1, 0, -1));
                float4 worldPos = mul(unity_CameraToWorld, viewPos);
                worldPos.xyz /= worldPos.w;

                float3 viewVector = normalize(_WorldSpaceCameraPos - worldPos.xyz);

                float3 p = _WorldSpaceCameraPos / kUnit - kEarthCenter;
                return float4(GetSkyRadiance(p, viewVector), 1.0);
            }
            ENDHLSL
        }
    }
}
