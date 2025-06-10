from mpi4py import MPI
import preprocessing
import algorithm_bootstrap
import algorithm
import time

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

start = time.time()

if rank == 0:
    preprocessing.preprocessing()
comm.Barrier()

if preprocessing.bootstrapping > 0:
    algorithm_bootstrap.algorithm()
else:
    algorithm.algorithm

if rank == 0:
    end = time.time()
    print(f"Runtime: {end - start} seconds")