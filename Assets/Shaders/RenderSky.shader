Shader "Hidden/Atmosphere"
{
    Properties
    {
    }

    SubShader
    {
        Tags { "RenderType"="Opaque" "RenderPipeline"="UniversalPipeline" }
        LOD 100
        Cull Off ZWrite Off ZTest Always

        Pass
        {
            HLSLPROGRAM
            #include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/Core.hlsl"
            #include "Packages/com.unity.render-pipelines.core/Runtime/Utilities/Blit.hlsl"

            #include "Definitions.hlsl"
            #include "UtilityFunctions.hlsl"
            #include "TransmittanceFunctions.hlsl"
            #include "ScatteringFunctions.hlsl"
            #include "RenderingFunctions.hlsl"
            
            #pragma vertex Vert
            #pragma fragment frag

            float3 earth_center;
            float3 sun_direction;

            float4 frag(Varyings input) : SV_Target
            {
                float4 viewPos = mul(unity_CameraInvProjection, float4(input.texcoord * 2 - 1, 0, -1));
                float4 worldPos = mul(unity_CameraToWorld, viewPos);
                worldPos.xyz /= worldPos.w;

                float3 viewVector = normalize(_WorldSpaceCameraPos - worldPos.xyz);

                float3 p = _WorldSpaceCameraPos - earth_center;

                return float4(GetSkyRadiance(p, viewVector, sun_direction), 1.0);
            }

            ENDHLSL
        }
    }
}
