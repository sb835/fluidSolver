using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class ParticleScript : MonoBehaviour
{
    private SimulationScript simulation;
    public int particle;
    public Vector2 particlePosition;
    public Vector2 particleVelocity;
    public Vector2 gridCell;
    public List<int> gridCellContent;
    public List<int> neighbors;
    public List<int> boundaryNeighbors;
    public float kernel;
    public Vector2 kernelGradient;
    public Vector2 kernelGradientSum;
    public int neighbor;
    public int oldPos;
    public float kernelSum;
    public float density;
    public float pressure;
    public Vector2 force;
    // Start is called before the first frame update
    void Start()
    {
        simulation = GameObject.FindGameObjectWithTag("Simulation").GetComponent<SimulationScript>();
        oldPos = neighbor;
    }

    // Update is called once per frame
    void Update()
    {
        particle = simulation.currentParticle;
        particlePosition = simulation.currentParticlePosition;
        particleVelocity = simulation.currentParticleVelocity;
        gridCell = simulation.currentGridCell;
        gridCellContent = simulation.contentCurrentGridCell;
        neighbors = simulation.neighbors[particle];
        boundaryNeighbors = simulation.boundaryNeighbors[particle];
        kernel = simulation.smoothingKernel(simulation.positions[particle], simulation.positions[neighbor], simulation.particleSpacing);
        kernelGradient = simulation.smoothingKernelDerivative(simulation.positions[particle], simulation.positions[neighbor], simulation.particleSpacing);
        // simulation.colors[neighbor] = Color.green;
        if (oldPos != neighbor)
        {
            simulation.colors[oldPos] = Color.blue;
            oldPos = neighbor;
        }
        float sum = 0.0f;
        Vector2 gradientSum = new Vector2(0, 0);
        foreach (int p in simulation.neighbors[particle])
        {
            sum += simulation.smoothingKernel(simulation.positions[particle], simulation.positions[p], simulation.particleSpacing);
            gradientSum += simulation.smoothingKernelDerivative(simulation.positions[particle], simulation.positions[p], simulation.particleSpacing);
        }
        kernelSum = sum;
        kernelGradientSum = gradientSum;
        density = simulation.densitys[particle];
        pressure = simulation.pressures[particle];
        force = simulation.forces[particle];
    }
}
