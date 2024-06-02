using System.Collections;
using System.Collections.Generic;
using Unity.VisualScripting;
using UnityEngine;

public class SimulationScript : MonoBehaviour
{
    public bool runSimulation = true;
    public Vector2[] positions;
    public Vector2[] velocitys;
    public Color[] colors;
    public List<int>[] neighbors;
    public int numParticles;
    public float particleSize;
    public float particleMass;
    public float particleSpacing;
    public float particleVolume;
    public float startDensity = 0.9f;
    public int texResolution = 2048;
    private float kernelSupportRadius;
    private Vector2 mouse;
    private DrawCirclesScript drawCirclesScript;
    private GridScript spatialGrid;
    public Camera mainCamera;
    private int previousParticle = 0;
    private Vector2 boundaries;


    // Forces
    public float gravity = -1;
    public float damping = 0.9f;


    // Not meant to be seen
    public int currentParticle;
    public Vector2 currentParticlePosition;
    public Vector2 currentParticleVelocity;
    public Vector2 currentGridCell;
    public List<int> contentCurrentGridCell;
    // Start is called before the first frame update
    void Start()
    {
        // Initialize all values
        drawCirclesScript = GameObject.FindGameObjectWithTag("DrawCircle").GetComponent<DrawCirclesScript>();
        spatialGrid = GameObject.FindGameObjectWithTag("Grid").GetComponent<GridScript>();
        texResolution = drawCirclesScript.texResolution;
        particleSpacing = spatialGrid.cellSize / 2;
        particleVolume = particleSpacing * particleSpacing;
        particleMass = startDensity * particleVolume;
        particleSize = particleMass;
        kernelSupportRadius = particleSpacing * 2;

        // Define boundaries of our simulation domain
        boundaries = new Vector2(19, 10);
        // InitializeParticles(32);
        // InitializeParticlesInCorners();
        spawnInEveryCell(3, 3, 25, 18);
        // fourParticlesInCell(new Vector2(10, 10), 0);
        drawCirclesScript.total = positions.Length;

        // Update grid
        spatialGrid.emptyGrid();
        spatialGrid.DrawGrid();

    }

    // Update is called once per frame
    void Update()
    {
        // Mouse input
        if (Input.GetMouseButtonDown(0))
        {
            ReturnMouse();
        }

        if (runSimulation)
        {
            SimulationStep(0.0166f);
        }

        DrawParticles(positions, colors);
    }


    void SimulationStep(float deltaTime)
    {
        // Update particle size
        particleMass = startDensity * particleVolume;
        particleSize = particleMass;

        // clear the grid
        spatialGrid.emptyGrid();

        // Add the number of each particle in the respective grid cell
        for (int i = 0; i < numParticles; i++)
        {
            Vector2 gridCoords = spatialGrid.computeCellPosition(positions[i]);
            if (spatialGrid.isValidCell(gridCoords))
            {
                spatialGrid.grid[(int)gridCoords.x, (int)gridCoords.y].Add(i);
            }
        }

        // Find all neighbors for each particle
        for (int i = 0; i < numParticles; i++)
        {
            findNeighbors(i);
        }
        // MoveParticles(deltaTime);
    }

    // Convert positions to texCoords
    private Vector2[] convertPositions(Vector2[] positions)
    {
        Vector2[] results = new Vector2[positions.Length];
        for (int i = 0; i < positions.Length; i++)
        {
            results[i] = positions[i] * 102.4f;
        }
        return results;
    }

    // Initialize one particle in each corner
    private void InitializeParticlesInCorners()
    {
        numParticles = 4;
        positions = new Vector2[4];
        colors = new Color[4];

        positions[0] = new Vector2(0, 0);
        colors[0] = Color.blue;

        positions[1] = new Vector2(0, boundaries.y);
        colors[1] = Color.blue;

        positions[2] = new Vector2(boundaries.x, 0);
        colors[2] = Color.blue;

        positions[3] = boundaries;
        colors[3] = Color.blue;
    }

    // Initialize particle
    private void InitParticle(int particleNum, Vector2 position, Vector2 velocity, Color color)
    {
        positions[particleNum] = position;
        velocitys[particleNum] = velocity;
        colors[particleNum] = color;
        neighbors[particleNum] = new List<int>();
    }
    // Initialize Particles
    private void InitializeParticles(int num)
    {
        numParticles = num;
        positions = new Vector2[num];
        velocitys = new Vector2[num];
        colors = new Color[num];
        neighbors = new List<int>[num];
        for (int i = 0; i < num; i++)
        {
            InitParticle(i, new Vector2(Random.Range(0, boundaries.x), Random.Range(0, boundaries.y)),
            new Vector2(0, 0),
            Color.blue);
        }
    }

    // Return mouse position
    private void ReturnMouse()
    {
        Vector3 screenPosition = Input.mousePosition;
        screenPosition.z = Camera.main.nearClipPlane + 1;
        mouse = Camera.main.ScreenToWorldPoint(screenPosition);

        for (int i = 0; i < numParticles; i++)
        {
            if (Vector2.Distance(positions[i], mouse) <= particleSize)
            {
                colors[previousParticle] = Color.blue;
                colorNeighbors(previousParticle, Color.blue);
                colorNeighbors(i, Color.yellow);
                colors[i] = Color.red;
                currentParticle = i;
                currentParticlePosition = positions[i];
                currentParticleVelocity = velocitys[i];
                currentGridCell = spatialGrid.computeCellPosition(currentParticlePosition);
                contentCurrentGridCell = spatialGrid.grid[(int)currentGridCell.x, (int)currentGridCell.y];
                previousParticle = i;
            }
        }
    }

    //Move particles downwards
    private void MoveParticles(float deltaTime)
    {
        for (int i = 0; i < positions.Length; i++)
        {
            positions[i] += velocitys[i] * deltaTime;
            velocitys[i] += Vector2.down * gravity * deltaTime;
            containParticles(i);
        }
    }

    // Draw all Particles at Position with their respective color
    private void DrawParticles(Vector2[] positions, Color[] colors)
    {
        drawCirclesScript.DrawCirclesAtPositions(convertPositions(positions), colors, particleSize * 102.4f);
        drawCirclesScript.DispatchKernel(Mathf.Max(numParticles / 16, 1));
    }

    // Neighbor search
    private void findNeighbors(int particle)
    {
        List<int> n = new List<int>();
        Vector2 gridCell = spatialGrid.computeCellPosition(positions[particle]);
        for (int x = -1; x <= 1; x++)
        {
            for (int y = -1; y <= 1; y++)
            {
                int cellX = (int)gridCell.x + x;
                int cellY = (int)gridCell.y + y;
                if (spatialGrid.isValidCell(new Vector2(cellX, cellY)))
                {
                    foreach (int p in spatialGrid.grid[cellX, cellY])
                    {
                        if (Vector2.Distance(positions[particle], positions[p]) < kernelSupportRadius)
                        {
                            n.Add(p);
                        }
                    }
                }
            }
        }
        neighbors[particle] = n;
    }

    // Testing 0: Helper functions

    private void containParticles(int particle)
    {
        if (positions[particle].x > boundaries.x)
        {
            positions[particle].x = boundaries.x;
            velocitys[particle] *= -1 * damping;
        }
        if (positions[particle].x < 0)
        {
            positions[particle].x = 0;
            velocitys[particle] *= -1 * damping;
        }
        if (positions[particle].y > boundaries.y)
        {
            positions[particle].y = boundaries.y;
            velocitys[particle] *= -1 * damping;
        }
        if (positions[particle].y < 0)
        {
            positions[particle].y = 0;
            velocitys[particle] *= -1 * damping;
        }
    }

    // My smoothing Kernel implemented as a cubic spline kernel
    public float smoothingKernel(Vector2 xi, Vector2 xj, float h)
    {
        float alpha = 5 / (14 * Mathf.PI * (h * h));
        float d = Vector2.Distance(xi, xj) / h;
        float t1 = Mathf.Max(1 - d, 0);
        float t2 = Mathf.Max(2 - d, 0);
        return alpha * (t2 * t2 * t2 - 4 * t1 * t1 * t1);
    }

    // My kernel Gradient
    Vector2 smoothingKernelDerivative(Vector2 xi, Vector2 xj, float h)
    {
        float alpha = 5 / (14 * Mathf.PI * (h * h));
        float d = Vector2.Distance(xi, xj) / h;
        float t1 = Mathf.Max(1 - d, 0);
        float t2 = Mathf.Max(2 - d, 0);
        if (Vector2.Distance(xi, xj) == 0)
        {
            return new Vector2(0, 0);
        }
        return alpha * (xi - xj) / (Vector2.Distance(xi, xj) * h) * (-3 * t2 * t2 + 12 * t1 * t1);
    }



    // Testing 1: Particle constants

    // Initiatlize 4 Particles in a grid cell
    void fourParticlesInCell(Vector2 gridCell, int particleNum)
    {
        // numParticles = 4;
        // positions = new Vector2[4];
        // velocitys = new Vector2[4];
        // colors = new Color[4];

        // Compute world coords of the upper right corner of our cell
        Vector2 coords = spatialGrid.computeWorldCoords((int)gridCell.x, (int)gridCell.y);
        float step = spatialGrid.cellSize / 4;
        Vector2 point1 = new Vector2(coords.x - step, coords.y - step);
        Vector2 point2 = new Vector2(coords.x - step, coords.y - 3 * step);
        Vector2 point3 = new Vector2(coords.x - 3 * step, coords.y - 3 * step);
        Vector2 point4 = new Vector2(coords.x - 3 * step, coords.y - step);
        InitParticle(particleNum, point1, new Vector2(0, 0), Color.blue);
        InitParticle(particleNum + 1, point2, new Vector2(0, 0), Color.blue);
        InitParticle(particleNum + 2, point3, new Vector2(0, 0), Color.blue);
        InitParticle(particleNum + 3, point4, new Vector2(0, 0), Color.blue);
    }

    // Spawn particles in every cell
    void spawnInEveryCell(int startX, int startY, int xDirection, int yDirection)
    {
        int num = xDirection * yDirection * 4;
        numParticles = num;
        positions = new Vector2[num];
        velocitys = new Vector2[num];
        colors = new Color[num];
        neighbors = new List<int>[num];
        int counter = 0;

        float boundX = Mathf.Min(spatialGrid.width - startX, startX + xDirection);
        float boundY = Mathf.Min(spatialGrid.height - startY, startY + yDirection);
        for (int x = startX + 1; x <= boundX; x++)
        {
            for (int y = startY + 1; y <= boundY; y++)
            {
                fourParticlesInCell(new Vector2(x, y), counter);
                counter += 4;
            }
        }
    }


    // Testing 2: Neighbor search

    // Color all neighbors in one color
    private void colorNeighbors(int particle, Color color)
    {
        foreach (int num in neighbors[particle])
        {
            colors[num] = color;
        }
    }
}
