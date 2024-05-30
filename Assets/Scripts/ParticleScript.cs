using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class ParticleScript : MonoBehaviour
{
    private SimulationScript simulation;
    public int particle;
    public Vector2 particlePositions;
    // Start is called before the first frame update
    void Start()
    {
        simulation = GameObject.FindGameObjectWithTag("Simulation").GetComponent<SimulationScript>();
    }

    // Update is called once per frame
    void Update()
    {
        particle = simulation.currentParticle;
        particlePositions = simulation.currentParticlePosition;
    }
}
