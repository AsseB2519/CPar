<div align="center">
  
# Parallel Computing

</div>

This project revolves around the optimization and analysis of a [Molecular Dynamics Simulation](https://github.com/FoleyLab/MolecularDynamics/tree/master) code designed for argon gas atoms. 

Our objective was to enhance the performance, particularly the execution time, of the simulation code through a series of optimization techniques and parallelization strategies.

<h3>Approaches:</h3>

- **Instruction-Level Parallelism (ILP) and Vectorization**: Optimizing the code to exploit parallel execution at the instruction level and leverage vector operations for processing multiple data elements simultaneously.
- **OpenMP**: Utilizing OpenMP directives to facilitate parallel execution in shared-memory systems, thereby distributing workload across multiple threads.
- **CUDA / MPI**: Leveraging CUDA for GPU acceleration and MPI for parallel programming across distributed-memory systems. This combination enables parallelism across multiple GPUs and/or nodes in a cluster.

<h3>Team Members:</h3>
<p> 
  
  - <a href="https://github.com/AsseB2519">Afonso Bessa</a>
  - <a href="https://github.com/CaraPokeBao">Francisco Claudino</a>
</p>

<div style="flex: 1;">
  <h3>Evaluation:</h3>
  <p><strong>Score:</strong> 17/20</p>
  <p><strong>Tools:</strong> C++ / ILP / Vectorization / OpenMP / MPI / CUDA</p>
</div>
