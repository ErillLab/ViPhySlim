from mpi4py import MPI
import preprocessing
import algorithm_bootstrap
import time

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

start = time.time()

if rank == 0:
    preprocessing.preprocessing()

if size > 1:   
    comm.Barrier()

algorithm_bootstrap.algorithm()

if rank == 0:
    end = time.time()
    print(f"Runtime: {end - start:.2f} seconds")