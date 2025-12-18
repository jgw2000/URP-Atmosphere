using UnityEngine;

public class RotateLight : MonoBehaviour
{
    public float speed = 5.0f;

    private Vector3 lastMousePos;

    private void Update()
    {
        Vector3 delta = (lastMousePos - Input.mousePosition) * Time.deltaTime * speed;

        if (Input.GetMouseButton(0) && Input.GetKey(KeyCode.LeftControl))
        {
            transform.Rotate(new Vector3(-delta.y, -delta.x, 0));
        }

        lastMousePos = Input.mousePosition;
    }
}
