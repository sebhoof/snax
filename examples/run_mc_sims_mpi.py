from mpi4py import MPI
import os
import sys
import time

script_path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, script_path+"/code")
from fluence_calc_mc import *

# Set up the MPI environment and variables.
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
ncores = comm.Get_size()

# Define output path and verbosity level.
output_path = "./"
verbosity_level = 0

# Define the grid of ALP masses log10(m/eV) and coupling log10(g/GeV^-1) and their combination.
logm = np.arange(3.0, 9.0, 0.05)
logg = np.arange(-13.0, -3.0, 0.05)
pairs = np.array(np.meshgrid(logm, logg)).T.reshape(-1,2)
ntasks = len(pairs)

start_time = time.time()

if (rank > 0):
    # Send first result to main process to signal that this proc is ready for more tasks
    lgm, lgg = pairs[rank-1]
    res = expected_photon_fluence(10**lgm,10**(lgg-9.0),verbosity_level)
    comm.Send(res, dest=0)
    for i in range(ntasks):
        task_id = comm.recv(source=0)
        if (task_id > ntasks):
            break
        lgm, lgg = pairs[task_id-1]
        res = expected_photon_fluence(10**lgm,10**(lgg-9.0),verbosity_level)
        comm.Send(res, dest=0)
    print('MPI rank {} finished! MC simulations took {:.1f} mins.'.format(rank, (time.time()-start_time)/60.0))

# Main process distribute tasks and receive results
if (rank == 0):
    print('Main process waiting for results from {} other processes...'.format(ncores-1), flush=True)
    all_results = []
    # Receive results and send more work (send out ncores work task too many to finalise jobs).
    for task_id in range(ncores, ntasks+ncores):
        info = MPI.Status()
        # The container for the result
        res = np.zeros(3)
        comm.Recv(res, source=MPI.ANY_SOURCE, status=info)
        all_results.append(res)
        worker_id = info.Get_source()
        comm.send(task_id, dest=worker_id)
    print('All MPI tasks finished after {:.1f} mins!'.format(rank, (time.time()-start_time)/60.0), flush=True)

    start_io_time = time.time()
    out_file_name = output_path+"SN1987A_DecayFluence.dat"
    print('Formatting results and saving them to '+out_file_name+'.')
    a = np.array(all_results)
    a = a[a[:,1].argsort()]
    a = a[a[:,0].argsort(kind='mergesort')]
    header_string = "Fluence of SN1987A photons from axion decay\nMC simulations by Marie Lecroq, Sebastian Hoof, and Csaba Balazs based on arXiv:1702.02964\nColumns: Value of log10(m/eV) | Value of log10(g/GeV^-1) | Fluence value in cm^-2"
    np.savetxt(out_file_name, a, fmt="%.3f %.3f %.10e", header=header_string, comments="# ")

    print('Formatting and saving file took {:.2f} mins!'.format((time.time()-start_io_time)/60.0), flush=True)
    print('All tasks complete! Finishing MPI routine now.', flush=True)
