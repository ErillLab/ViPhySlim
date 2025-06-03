from mpi4py import MPI
import preprocessing
import algorithm

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

if rank == 0:
    preprocessing.preprocessing()
comm.Barrier()

algorithm.algorithm()