import csv
import math
import matplotlib.pyplot as plt
import multiprocessing
import numpy as np
import os
import sys
import time
import uproot
import click
from tqdm import tqdm

def partition_helper(slice_entries, file_entries, file_curr, entry_curr):
    if slice_entries <= file_entries[file_curr] - entry_curr:
        return [file_curr, slice_entries + entry_curr]
    elif file_curr == len(file_entries) - 1:
        return [file_curr, file_entries[-1]]
    else:
        return partition_helper(slice_entries - file_entries[file_curr] + entry_curr, file_entries, file_curr + 1, 0)

def partition(files, n_processes):
    file_entries = [file.num_entries for file in files]
    slice_entries = math.ceil(sum(file_entries) / n_processes)
    slices = []
    file_start = 0
    entry_start = 0
    for i in range(n_processes):
        slices.append([file_start, entry_start] + partition_helper(slice_entries, file_entries, file_start, entry_start))
        file_start = slices[-1][-2]
        entry_start = slices[-1][-1]
    return slices

def read_slice(files, slices, index, data):
    data_slice = []
    for i in range(slices[index][0], slices[index][2] + 1):
        data_slice.append(files[i].arrays("candidate_vMass",
                              "(candidate_charge == 0)\
                              & (candidate_cosAlpha > 0.99)\
                              & (candidate_lxy / candidate_lxyErr > 3.0)\
                              & (candidate_vProb > 0.05)\
                              & (ditrack_mass > 1.014) & (ditrack_mass < 1.024)\
                              & (candidate_vMass > 5.33) & (candidate_vMass < 5.4)",
                              entry_start=slices[index][1] if i == slices[index][0] else None,
                              entry_stop=slices[index][3] if i == slices[index][2] else None,
                              array_cache=None,
                              library="np")["candidate_vMass"])
    data.append(np.concatenate(tuple(data_slice)))

def runtime_measure_mt(path, n_files, n_processes):
    if n_files == 0: return 0
    if n_processes == 0: return runtime_measure(path, n_files)
    start = time.time()
    files = [uproot.open(path=path + filename + ":rootuple/CandidateTree", object_cache=None, array_cache=None) for filename in sorted(os.listdir(path))[:n_files]]
    slices = partition(files, n_processes)
    data = multiprocessing.Manager().list()
    processes = []
    for i in range(n_processes):
        p = multiprocessing.Process(target=read_slice, args=[files, slices, i, data])
        p.start()
        processes.append(p)

    for p in processes:
        p.join()

    np.concatenate(tuple(data))

    return time.time() - start

def runtime_measure(path, n_files):
    if n_files == 0: return 0
    start = time.time()
    files = [uproot.open(path=path + filename + ":rootuple/CandidateTree", object_cache=None, array_cache=None) for filename in sorted(os.listdir(path))[:n_files]]
    data = []
    for file in files:
        data.append(file.arrays("candidate_vMass",
                              "(candidate_charge == 0)\
                              & (candidate_cosAlpha > 0.99)\
                              & (candidate_lxy / candidate_lxyErr > 3.0)\
                              & (candidate_vProb > 0.05)\
                              & (ditrack_mass > 1.014) & (ditrack_mass < 1.024)\
                              & (candidate_vMass > 5.33) & (candidate_vMass < 5.4)",
                              array_cache=None,
                              library="np")["candidate_vMass"])

    np.concatenate(tuple(data))

    return time.time() - start

def runtime_vs_variable(path, target_dir, measure_function, variable, step, n_loops, var_max, constant=None):
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)

    result_path = ("%s/runtime_vs_%s_%d_%d_%d_%d.csv" % (target_dir, variable, constant, var_max, step, n_loops)) if constant else ("%s/runtime_vs_%s_%d_%d_%d.csv" % (target_dir, variable, var_max, step, n_loops))

    x = [a for a in range(0, var_max + step, step)]

    if not os.path.exists(result_path):
        with open(result_path, "w", newline="") as f:
            csv.writer(f).writerow(x)
    print(">Warming up. First measurement will be ignored.")
    y = [measure_function(*(path, i if "size" in variable else constant, constant if "size" in variable else i) if constant else (path, i)) for i in range(0, step, step)]

    for n in tqdm(range(n_loops)):
        with open(result_path, "r") as f:
            if sum(1 for row in csv.reader(f)) == n_loops + 1: break
        y = [measure_function(*(path, i if "size" in variable else constant, constant if "size" in variable else i) if constant else (path, i)) for i in range(0, var_max + step, step)]
        with open(result_path, "a", newline="") as f:
            csv.writer(f).writerow(y)

usefullvars = ["size","size_mt","threads"]

@click.command()
@click.option('--step', default=4, help='Scanning step.')
@click.option('--loops', default=20, help='No. of measurements to do.')
@click.option('--max', default=128, help='Variable max.')
@click.option('--name', prompt='Path', default = "../data/32_files/",
              help='File path.')
@click.option('--output', prompt='Output', default = "runtime_tests_uproot/32_files/",
              help='Output dir.')
@click.option('--variable', prompt='Variable to test', default = "processes",
              help='Running variable (either size_mp, size or threads).')
@click.option('--n_threads', prompt='Number of threads', default = 32,
              help='Number of threads to run. If variable==threads this will be useless.')
@click.option('--n_files', prompt='Number of files', default = 128,
              help='Number of files to load (used when testing multiple processes variable).')
def run(step,loops,max,name,output,variable,n_threads,n_files):

    if variable not in usefullvars:
        raise Exception("Sorry, only %s as variables. No %s here!"%(usefullvars,variable))

    if variable == "size":
        runtime_vs_variable(path, target_dir, runtime_measure, variable, step, loops, max)
    elif variable == "size_mt":
        runtime_vs_variable(path, target_dir, runtime_measure_mt, "sivariable, step, loops, max, n_threads)
    elif variable == "threads":
        runtime_vs_variable(path, target_dir, runtime_measure_mt, variable, step, loops, max, n_files)

if __name__ == '__main__':
    run()
