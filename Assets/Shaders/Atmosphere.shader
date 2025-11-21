Shader "Custom/Atmosphere"
{
    Properties
    {
        [IntRange] _InScatteringPoints("In Scattering Points", Range(5, 15)) = 10
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
                // this is needed so we account XR platform differences in how they handle texture arrays
                UNITY_SETUP_STEREO_EYE_INDEX_POST_VERTEX(input);

                float4 viewPos = mul(unity_CameraInvProjection, float4(input.texcoord * 2 - 1, 0, -1));
                float4 worldPos = mul(unity_CameraToWorld, viewPos);
                worldPos.xyz /= worldPos.w;

                float3 viewVector = normalize(_WorldSpaceCameraPos - worldPos.xyz);

                if (viewVector.y > 0)
                    return float4(0.5, 0.7, 1.0, 1);
                else
                    return float4(CalculateAtmosphere(float3(0, Rg, 0), viewVector), 1);
            }
            ENDHLSL
        }
    }
}
