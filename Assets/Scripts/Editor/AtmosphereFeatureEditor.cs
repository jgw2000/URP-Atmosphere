using UnityEngine;
using UnityEditor;

[CustomEditor(typeof(AtmosphereFeature))]
public class AtmosphereFeatureEditor : Editor
{
    public override void OnInspectorGUI()
    {
        base.OnInspectorGUI();

        if (GUILayout.Button("Bake Transmittance LUT"))
        {
            AtmosphereFeature feature = (AtmosphereFeature)target;
            feature.BakeTransmittanceLUT();
        }
    }
}
