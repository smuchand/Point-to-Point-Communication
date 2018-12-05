# Point-to-Point-Communication
This repository is based on a course project which focuses on MPI Point-to-Point Communication. 
There are 3 examples in this repository:

1) Numerical Integration Problem using Master-Worker model:

      A master-worker system enables dynamic scheduling of applications. The idea is that one of the MPI processes (usually rank 0) will be responsible for giving work to the other MPI processes and collect results. Usually the master process starts by sending one chunk of work to all the workers and when one of the worker provides the result of that chunk of work, the master process will send a new chunk of work to the worker that just completed the work. The workers will wait for a chunk of work to perform, perform that chunk of work, and return the result to the master node. Also, the master node needs to notify the worker nodes when there is no more work to perform so the workers can quit gracefully. 
  
  
2) Heat 2D:

      The problem is defined on a discrete 2D space of size n × n; let’s call it H. Initialize H<sub>0</sub> using the provided functions. The kth iteration of the heat equation is defined by H<sub>k</sub> is defined by:
          
        
      H<sub>k</sub> [i][j] = 1/5[(H<sub>k-1</sub> [i][j]) + (H<sub>k-1</sub> [i-1][j]) + (H<sub>k-1</sub> [i][j-1]) + (H<sub>k-1</sub> [i][j+1]) +  (H<sub>k-1</sub> [i+1][j])] 


3) Numerical Integration using Advanced Master-Worker model:

    One of the issue with simple master-worker implementation is that the worker process will be idle until the master provides the next chunk. To avoid this, it is common that the master will give more than a single chunk of work (say, 3 chunks of work) to each worker. That way, when a worker sends a result to the master node, the worker still has some chunks of work (here, 2 chunks) to perform.
