using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class SimulationScript : MonoBehaviour
{
    public Vector2[] positions;
    public Color[] colors;
    public int numParticles;
    public float particleSize;
    private int texResolution;
    private Vector2 mouse;
    private DrawCirclesScript drawCirclesScript;
    public Camera mainCamera;
    private int previousParticle = 0;


    // Not meant to be seen
    public int currentParticle;
    public Vector2 currentParticlePosition;
    // Start is called before the first frame update
    void Start()
    {
        // Initialize all values
        drawCirclesScript = GameObject.FindGameObjectWithTag("DrawCircle").GetComponent<DrawCirclesScript>();
        texResolution = drawCirclesScript.texResolution;
        particleSize = 50.0f;
        InitializeParticle(32);
        drawCirclesScript.total = positions.Length;

        // Define boundaries of our simulation domain

    }

    // Update is called once per frame
    void Update()
    {
        DrawParticles(positions, colors, particleSize);
        // MoveParticle();
        if (Input.GetMouseButtonDown(0))
        {
            returnMouse();
        }
    }

    // Initialize Particle
    private void InitializeParticle(int num)
    {
        numParticles = num;
        positions = new Vector2[num];
        colors = new Color[num];
        for (int i = 0; i < num; i++)
        {
            positions[i] = new Vector2(Random.Range(0, 2048), Random.Range(0, 2048));
            colors[i] = Color.blue;
        }
    }

    // Return mouse position
    private void returnMouse()
    {
        RaycastHit hit;
        Ray ray = mainCamera.ScreenPointToRay(Input.mousePosition);

        if (Physics.Raycast(ray, out hit))
        {
            mouse.x = hit.textureCoord.x * texResolution;
            mouse.y = hit.textureCoord.y * texResolution;
        }
        for (int i = 0; i < numParticles; i++)
        {
            if (Vector2.Distance(positions[i], mouse) <= particleSize)
            {
                colors[previousParticle] = Color.blue;
                colors[i] = Color.red;
                currentParticle = i;
                currentParticlePosition = positions[i];
                previousParticle = i;
            }
        }
    }

    //Move particles downwards
    private void MoveParticle()
    {
        for (int i = 0; i < positions.Length; i++)
        {
            positions[i].y -= 1.0f;
        }
    }

    // Draw all Particles at Position with their respective color
    private void DrawParticles(Vector2[] positions, Color[] colors, float radius)
    {
        drawCirclesScript.DrawCirclesAtPositions(positions, colors, particleSize);
        drawCirclesScript.DispatchKernel(numParticles / 32);
    }
}
