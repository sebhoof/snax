from mpi4py import MPI

import time
import numpy as np
from datetime import date


def task_wrapper(task, id, pairs, ebins, tbins):
    ie, lgm, lgg = pairs[id]
    ebin = ebins[int(ie)]
    t0 = time.time()
    res = task(10**lgm, 10**lgg, ebin, tbins)
    t1 = time.time()
    dt = (t1-t0)/60.0
    res = np.array(list(ebin) + [lgm, lgg] + res + [dt])
    return res, dt

start_time = time.time()

def run_mpi_job(task, out_file_name_root, ebins, tbins, lgmvals, lggvals, save_temp_results=True, save_timing_info=True):
    # Set up the MPI environment and variables.
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    ncores = comm.Get_size()

    ntbins = len(tbins)-1
    nebins = len(ebins)
    pairs = np.array(np.meshgrid(range(nebins), lgmvals, lggvals)).T.reshape(-1,3)
    ntasks = len(pairs)
    
    if (rank > 0):
        # Send first result to main process to signal that this proc is ready for more tasks
        res, dt = task_wrapper(task, rank-1, pairs, ebins, tbins)
        print('Initial job finished on rank {}, took {:.2} mins.'.format(rank, dt), flush=True)
        comm.Send(res, dest=0)
        for i in range(ntasks):
            task_id = comm.recv(source=0)
            if (task_id > ntasks):
                break
            res, dt = task_wrapper(task, task_id-1, pairs, ebins, tbins)
            comm.Send(res, dest=0)
        dt = (time.time()-start_time)/3600.0
        print('MPI rank {} finished! MC simulations took {:.3f} hrs.'.format(rank, dt))
        
    # Main process distribute tasks and receive results
    if (rank == 0):
        all_results = []
        out_file_name = out_file_name_root+".dat"
        out_file_temp_name = out_file_name_root+"_temp.dat"
        n_info = 1
        buffer_size = 4+ntbins+1
        n_temp_save = int(ntasks//20 + 1) # Save after ~5% progress for many tasks (otherwise every task)
        print('Main process waiting for {} results from {} other processes...'.format(ntasks, ncores-1), flush=True)

        # Receive results and send more work (send out ncores work task too many to finalise jobs).
        for task_id in range(ncores, ntasks+ncores):
            info = MPI.Status()
            res = np.zeros(buffer_size) # Container for the result
            comm.Recv(res, source=MPI.ANY_SOURCE, status=info)
            worker_id = info.Get_source()
            if save_timing_info:
                all_results.append(list(res))
            else:
                all_results.append(list(res[:-1]))
            if save_temp_results and (task_id%n_temp_save == 0):
                t0 = time.time()
                a = np.array(all_results)
                np.savetxt(out_file_temp_name, a, fmt="%.6e")
                t1 = time.time()
                dt = t1 - t0
                print('Job {}/{} finished on rank {}. Writing temporary results took {:.1f} s.'.format(task_id-ncores+1, ntasks, worker_id, dt), flush=True)
            comm.send(task_id, dest=worker_id)
        print('All MPI tasks finished after {:.3f} hrs!'.format( (time.time()-start_time)/3600.0 ), flush=True)
        print('Formatting results and saving them to '+out_file_name+'...')
        a = np.array(all_results)
        # Sort by energy bin (lower end, col 0), then mass (col 2), then coupling (col 3)
        a = a[np.lexsort((a[:,3],a[:,2],a[:,0]))]
        header = "ALP fluence from supernova SN1987A, calculated on {}\n".format(date.today().strftime("%Y-%m-%d"))
        header += "Energy bins: {}\n".format(ebins)
        header += "Time bins: {}\n".format(tbins)
        check = (ntbins > 1)
        header += "Columns: (1-2) Energy bin [eV] | (3) Mass log10(m/eV) | (4) Coupling log10(g/eV^-1) | (5{}) Fluence in time bin{} [cm^-2]".format("-{}".format(buffer_size-1) if check else "", "s" if check else "")
        if save_timing_info:
            header += " | ({}) Computation time [min]".format(buffer_size)
        np.savetxt(out_file_name, a, fmt="%.6e", header=header)
        print('All tasks complete! Finishing MPI routine...', flush=True)