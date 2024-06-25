using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using UnityEngine;

public class SimulationScript : MonoBehaviour
{
    public bool moveParticles = false;
    public bool colorParticles = false;
    public bool withSplitting = true;
    public bool countAvgDensity = false;
    public int tests = 1;
    public List<float> averageDensity;
    public Vector2 amountParticles;
    public Vector2 particlePosition;
    public float particleStartSpacing;
    public float speedColor = 4.0f;
    public float alpha;
    public Vector2[] positions;
    public Vector2[] velocitys;
    public Color[] colors;
    public List<int>[] neighbors;
    public List<int>[] boundaryNeighbors;
    public float[] densitys;
    public float[] pressures;
    public Vector2[] forces;
    public Vector2[] nPForces;
    public Vector2[] predictedPositions;
    public Vector2[] predictedVelocitys;
    public Vector2[] boundaryPositions;
    public Color[] boundaryColors;
    public Vector2 maxVelocity = new Vector2(2, 9);
    public float timeStepMultiplyer = 0.9f;
    public float timeStep;
    public float stiffness = 100.0f;
    public float v;
    public int numParticles;
    public int numBoundaries;
    public float particleSize;
    public float particleMass;
    public float particleSpacing;
    public float particleVolume;
    public float startDensity = 1.5f;
    public int texResolution = 2048;
    public float kernelSupportRadius;
    private DrawCirclesScript drawCirclesScript;
    private GridScript spatialGrid;
    public Camera mainCamera;
    private Vector2 mouse;
    private Vector2 start;
    private Vector2 boundaries;


    // Forces
    public float gravity = -1;


    // Not meant to be seen
    public int currentParticle;
    public Vector2 currentGridCell;
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
        start = new Vector2(1, 1);
        boundaries = new Vector2(16, 9);

        initializeParticles((int)amountParticles.x, (int)amountParticles.y, particlePosition, particleStartSpacing);

        // Update grid
        spatialGrid.emptyGrid();
        spatialGrid.DrawGrid();

        // Initialize boundary particles
        initializeBorder();

        numParticles = positions.Length;
        numBoundaries = boundaryPositions.Length;

        // Inform the shader about the total amount of drawn particles
        drawCirclesScript.total = positions.Length + boundaryPositions.Length;

        averageDensity = new List<float>();

    }
    void ResetValues(int amount)
    {
        positions = new Vector2[amount];
        velocitys = new Vector2[amount];
        colors = new Color[amount];
        densitys = new float[amount];
        pressures = new float[amount];
        forces = new Vector2[amount];
        nPForces = new Vector2[amount];
        predictedPositions = new Vector2[amount];
        predictedVelocitys = new Vector2[amount];
        boundaryPositions = new Vector2[amount];
        boundaryColors = new Color[amount];
        neighbors = new List<int>[amount];
        boundaryNeighbors = new List<int>[amount];
    }

    void initializeParticles(int numX, int numY, Vector2 start, float spacing)
    {
        int amount = numX * numY;
        int index = 0;
        ResetValues(amount);
        for (int x = 0; x < numX; x++)
        {
            float xx = start.x + spacing * x;
            for (int y = 0; y < numY; y++)
            {
                float yy = start.y + spacing * y;
                InitParticle(index, new Vector2(xx, yy), Color.blue);
                index++;
            }
        }

    }
    void initializeBorder()
    {
        int index = 0;
        // Left row
        for (float y = start.y; y < boundaries.y; y += particleSize * 2)
        {
            boundaryPositions[index] = new Vector2(start.x, y);
            boundaryColors[index] = Color.black;
            index++;
        }


        // Double
        for (float y = start.y; y < boundaries.y + 4 * particleSize; y += particleSize * 2)
        {
            boundaryPositions[index] = new Vector2(start.x - 2 * particleSize, y);
            boundaryColors[index] = Color.black;
            index++;
        }

        // Right row
        for (float y = start.y; y < boundaries.y; y += particleSize * 2)
        {
            boundaryPositions[index] = new Vector2(boundaries.x, y);
            boundaryColors[index] = Color.black;
            index++;
        }

        // Double
        for (float y = start.y - 2 * particleSize; y < boundaries.y + 2 * particleSize; y += particleSize * 2)
        {
            boundaryPositions[index] = new Vector2(boundaries.x + 2 * particleSize, y);
            boundaryColors[index] = Color.black;
            index++;
        }

        // Bottom row
        for (float x = start.x; x < boundaries.x; x += particleSize * 2)
        {
            boundaryPositions[index] = new Vector2(x, start.y);
            boundaryColors[index] = Color.black;
            index++;
        }

        // Double
        for (float x = start.x - 2 * particleSize; x < boundaries.x; x += particleSize * 2)
        {
            boundaryPositions[index] = new Vector2(x, start.y - 2 * particleSize);
            boundaryColors[index] = Color.black;
            index++;
        }

        // Top row
        for (float x = start.x; x < boundaries.x; x += particleSize * 2)
        {
            boundaryPositions[index] = new Vector2(x, boundaries.y + 0.5f * particleSize);
            boundaryColors[index] = Color.black;
            index++;
        }

        // Double
        for (float x = start.x; x < boundaries.x + 2 * particleSize; x += particleSize * 2)
        {
            boundaryPositions[index] = new Vector2(x, boundaries.y + 2.5f * particleSize);
            boundaryColors[index] = Color.black;
            index++;
        }
    }

    // Update is called once per frame
    void Update()

    {
        // Start Stop Simulation
        if (Input.GetKeyDown(KeyCode.Space))
        {
            moveParticles = !moveParticles;
        }

        // Restart
        if (Input.GetKeyDown(KeyCode.Backspace))
        {
            Start();
        }

        // Clean particles
        if (Input.GetKeyDown(KeyCode.C))
        {
            ResetValues(0);
            numParticles = 0;
            numBoundaries = 0;
        }

        if (Input.GetKey(KeyCode.S))
        {
            Vector3 screenPosition = Input.mousePosition;
            screenPosition.z = Camera.main.nearClipPlane + 1;
            mouseForceOut(Camera.main.ScreenToWorldPoint(screenPosition), 2f, 1.5f);
        }

        if (Input.GetKey(KeyCode.D))
        {
            Vector3 screenPosition = Input.mousePosition;
            screenPosition.z = Camera.main.nearClipPlane + 1;
            mouseForceIn(Camera.main.ScreenToWorldPoint(screenPosition), 2f, 1.5f);
        }

        // Mouse input
        if (Input.GetMouseButton(0))
        {
            ReturnMouse();
        }
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

        // Add the number of each boundary particle in a gird cell
        // We distinguish the boundary particle from regular particles
        // by adding numParticles to their particle index
        for (int i = 0; i < numBoundaries; i++)
        {
            Vector2 gridCoords = spatialGrid.computeCellPosition(boundaryPositions[i]);
            if (spatialGrid.isValidCell(gridCoords))
            {
                spatialGrid.grid[(int)gridCoords.x, (int)gridCoords.y].Add(i + numParticles);
            }
        }
        // Check for max velocity
        Parallel.For(0, numParticles, i =>
        {
            if (velocitys[i].magnitude > maxVelocity.magnitude && velocitys[i].y > -15 && velocitys[i].magnitude > 0)
            {
                maxVelocity = velocitys[i];
            }
        });

        if (countAvgDensity)
        {
            // Add average density
            float avgDensity = 0.0f;
            Parallel.For(0, numParticles, i =>
            {
                avgDensity += densitys[i];
            });
            averageDensity.Add(avgDensity / numParticles);
        }

        // Change color according to velocity
        if (colorParticles)
        {
            Parallel.For(0, numParticles, i =>
            {
                float speed = velocitys[i].magnitude / speedColor;
                // float speed = 1.15f;
                colors[i] = Color.Lerp(Color.blue, Color.red, speed);
            });
        }
        // Calculate time step
        // timeStep = timeStepMultiplyer * (particleSpacing / maxVelocity.magnitude);
        SimulationStep(0.02f);

        DrawParticles();

    }
    void SimulationStep(float deltaTime)
    {

        // Find all neighbors for each particle
        Parallel.For(0, numParticles, i =>
        {
            findNeighbors(i);
        });

        // Compute non-pressure accelerations
        Parallel.For(0, numParticles, i =>
        {
            Vector2 v = computeViscosityAcceleration(i);
            Vector2 g = new Vector2(0, gravity);
            nPForces[i] = v + g;
        });

        // Predict next positions
        Parallel.For(0, numParticles, i =>
        {
            if (withSplitting)
            {
                predictedVelocitys[i] = velocitys[i] + nPForces[i] * deltaTime;
                predictedPositions[i] = positions[i] + predictedVelocitys[i] * deltaTime;
            }
            else
            {
                predictedVelocitys[i] = velocitys[i];
                predictedPositions[i] = positions[i];
            }
        });

        // Compute densitys
        Parallel.For(0, numParticles, i =>
        {
            computeDensity(i);
        });
        // Compute pressures
        Parallel.For(0, numParticles, i =>
        {
            pressures[i] = Mathf.Max(stiffness * ((densitys[i] / startDensity) - 1), 0);
        });

        // Compute pressure forces
        Parallel.For(0, numParticles, i =>
        {
            computePressureAcceleration(i);
        });

        // Update particle
        Parallel.For(0, numParticles, i =>
        {
            Vector2 acceleration = new Vector2(0, 0);
            if (withSplitting)
            {
                acceleration = forces[i];
            }
            else
            {
                acceleration = nPForces[i] + forces[i];
            }
            if (moveParticles)
            {
                MoveParticles(acceleration, i, deltaTime);
            }
        });
    }

    //Move particles
    private void MoveParticles(Vector2 acceleration, int particle, float deltaTime)
    {
        velocitys[particle] = predictedVelocitys[particle] + acceleration * deltaTime;
        // velocitys[particle] += acceleration * deltaTime;
        positions[particle] = positions[particle] + velocitys[particle] * deltaTime;
    }

    // Draw all Particles at Position with their respective color
    private void DrawParticles()
    {
        Vector2[] particlePositions = positions.Concat(boundaryPositions).ToArray();

        Color[] particleColors = colors.Concat(boundaryColors).ToArray();

        drawCirclesScript.DrawCirclesAtPositions(convertPositions(particlePositions), particleColors, particleSize * 102.4f);
        drawCirclesScript.DispatchKernel(Mathf.Max(particlePositions.Length / 16, 1));
    }

    // Convert positions to texCoords
    private Vector2[] convertPositions(Vector2[] pos)
    {
        Vector2[] results = new Vector2[pos.Length];
        for (int i = 0; i < pos.Length; i++)
        {
            results[i] = pos[i] * 102.4f;
        }
        return results;
    }

    // Initialize particle
    private void InitParticle(int index, Vector2 position, Color color)
    {
        positions[index] = position;
        velocitys[index] = new Vector2(0, 0);
        colors[index] = color;
        densitys[index] = 0.0f;
        pressures[index] = 0.0f;
        forces[index] = new Vector2(0, 0);
        nPForces[index] = new Vector2(0, 0);
        predictedPositions[index] = new Vector2(0, 0);
        predictedVelocitys[index] = new Vector2(0, 0);
    }
    // Initialize Particles
    private void InitializeRandomParticles(int num)
    {
        int index = 0;
        for (int i = 0; i < num; i++)
        {
            InitParticle(index, new Vector2(Random.Range(start.x + 1, boundaries.x - 1), Random.Range(start.y + 1, boundaries.y - 1)),
            Color.blue);
        }
    }

    // Initialize square of particles
    private void spawnSquareParticles(Vector2 origin, int width, int height, float particleSpace)
    {
        int index = 0;
        for (float x = origin.x; x < origin.x + (width) * 2 * (particleSize + particleSpace); x += 2 * (particleSize + particleSpace))
        {
            for (float y = origin.y; y < origin.y + (height - 1) * 2 * (particleSize + particleSpace); y += 2 * (particleSize + particleSpace))
            {
                InitParticle(index, new Vector2(x, y), Color.blue);
            }
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
                currentParticle = i;
                currentGridCell = spatialGrid.computeCellPosition(positions[i]);
            }
        }
    }

    // Neighbor search
    private void findNeighbors(int particle)
    {
        List<int> n = new List<int>();
        List<int> b = new List<int>();
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
                        if (p >= numParticles)
                        {
                            if (Vector2.Distance(positions[particle], boundaryPositions[p - numParticles]) < kernelSupportRadius)
                            {
                                b.Add(p - numParticles);
                            }
                        }
                        else
                        {
                            if (Vector2.Distance(positions[particle], positions[p]) < kernelSupportRadius)
                            {
                                n.Add(p);
                            }
                        }
                    }
                }
            }
        }
        neighbors[particle] = n;
        boundaryNeighbors[particle] = b;
    }


    // Compute viscosity acceleration of particle i
    private Vector2 computeViscosityAcceleration(int particle)
    {
        Vector2 viscosity = new Vector2(0, 0);
        foreach (int num in neighbors[particle])
        {
            Vector2 xij = positions[particle] - positions[num];
            Vector2 vij = velocitys[particle] - velocitys[num];
            Vector2 formula = vij * xij / ((xij * xij) + new Vector2(0.01f * particleSpacing * particleSpacing, 0.01f * particleSpacing * particleSpacing));
            Vector2 gradient = smoothingKernelDerivative(positions[particle], positions[num], particleSpacing);
            viscosity += particleMass / densitys[num] * formula * gradient;
        }
        return 2 * v * viscosity;
    }

    // Compute density for each particle
    private void computeDensity(int particle)
    {
        float result = 0.0f;
        foreach (int num in neighbors[particle])
        {
            result += particleMass * smoothingKernel(predictedPositions[particle], predictedPositions[num], particleSpacing);
        }

        foreach (int num in boundaryNeighbors[particle])
        {
            result += particleMass * smoothingKernel(predictedPositions[particle], boundaryPositions[num], particleSpacing);
        }
        densitys[particle] = result;
    }

    // Compute pressure acceleration
    private void computePressureAcceleration(int particle)
    {
        Vector2 result = new Vector2(0, 0);
        foreach (int num in neighbors[particle])
        {
            Vector2 gradient = smoothingKernelDerivative(positions[particle], positions[num], particleSpacing);
            float formula = (pressures[particle] / (densitys[particle] * densitys[particle])) + (pressures[num] / (densitys[num] * densitys[num]));
            result += particleMass * formula * gradient;
        }

        Vector2 result2 = new Vector2(0, 0);
        foreach (int num in boundaryNeighbors[particle])
        {
            Vector2 gradient = smoothingKernelDerivative(positions[particle], boundaryPositions[num], particleSpacing);
            float formula = (pressures[particle] / (densitys[particle] * densitys[particle])) + (pressures[particle] / (densitys[particle] * densitys[particle]));
            result2 += particleMass * formula * gradient;
        }
        forces[particle] = -result - result2;
    }

    // My smoothing Kernel implemented as a cubic spline kernel
    public float smoothingKernel(Vector2 xi, Vector2 xj, float h)
    {
        float d = Vector2.Distance(xi, xj) / h;
        float t1 = Mathf.Max(1 - d, 0);
        float t2 = Mathf.Max(2 - d, 0);
        return alpha * (t2 * t2 * t2 - 4 * t1 * t1 * t1);
    }

    // My kernel Gradient
    public Vector2 smoothingKernelDerivative(Vector2 xi, Vector2 xj, float h)
    {
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
    void fourParticlesInCell(Vector2 gridCell)
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
        InitParticle(0, point1, Color.blue);
        InitParticle(0, point2, Color.blue);
        InitParticle(0, point3, Color.blue);
        InitParticle(0, point4, Color.blue);
    }

    // Spawn particles in every cell
    void spawnInEveryCell(int startX, int startY, int xDirection, int yDirection)
    {
        float boundX = Mathf.Min(spatialGrid.width - startX, startX + xDirection);
        float boundY = Mathf.Min(spatialGrid.height - startY, startY + yDirection);
        for (int x = startX + 1; x <= boundX; x++)
        {
            for (int y = startY + 1; y <= boundY; y++)
            {
                fourParticlesInCell(new Vector2(x, y));
            }
        }
    }


    // Testing 2: Neighbor search

    // Color all neighbors in one color
    public void colorNeighbors(int particle, Color color)
    {
        foreach (int num in neighbors[particle])
        {
            colors[num] = color;
        }
    }

    public void colorBoundaryNeighbors(int particle, Color color)
    {
        foreach (int num in boundaryNeighbors[particle])
        {
            boundaryColors[num] = color;
        }
    }

    // Testing Mouse Controls

    // Push particles away from the cursor
    void mouseForceOut(Vector2 mousePos, float radius, float strength)
    {
        for (int i = 0; i < numParticles; i++)
        {
            if (Vector2.Distance(mousePos, positions[i]) < radius)
            {
                Vector2 direction = mousePos - positions[i];
                direction = direction.normalized;
                velocitys[i] = -direction * strength;
            }
        }
    }

    // Pull particles to the cursor
    void mouseForceIn(Vector2 mousePos, float radius, float strength)
    {
        for (int i = 0; i < numParticles; i++)
        {
            if (Vector2.Distance(mousePos, positions[i]) < radius)
            {
                Vector2 direction = mousePos - positions[i];
                direction = direction.normalized;
                velocitys[i] = direction * strength;
            }
        }
    }
}
