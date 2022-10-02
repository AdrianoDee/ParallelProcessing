import csv
import matplotlib.pyplot as plt
import numpy as np
import os
import ROOT
import sys
import time
import uproot
import click
from tqdm import tqdm


def partition(path, n_files, n_processes):
    filenames = sorted(os.listdir(path))
    partitions = []
    curr = 0
    for i in range(n_processes):
        if i >= n_files: break
        n_files_in_partition = n_files // n_processes if i >= n_files % n_processes else n_files // n_processes + 1
        files_to_read = ROOT.std.vector('string')()
        for j in range(n_files_in_partition):
            files_to_read.push_back(path + filenames[curr + j])
        curr += n_files_in_partition
        partitions.append(files_to_read)
    return partitions

def to_numpy(files, result):
    if files.empty(): return
    data = ROOT.RDataFrame("rootuple/CandidateTree", files)
    cut = data.Filter("candidate_charge == 0")\
          .Filter("candidate_cosAlpha > 0.99")\
          .Filter("candidate_vProb > 0.05")\
          .Filter("candidate_lxy / candidate_lxyErr > 3.0")\
          .Filter("ditrack_mass > 1.014")\
          .Filter("ditrack_mass < 1.024")\
          .Filter("candidate_vMass > 5.33")\
          .Filter("candidate_vMass < 5.40")
    result.append(cut.AsNumpy(["candidate_vMass"])["candidate_vMass"])

def runtime_measure(path, n_files, mt):

    if n_files == 0: return 0
    if mt:
        ROOT.ROOT.EnableImplicitMT()
    else:
        ROOT.ROOT.DisableImplicitMT()

    filenames = sorted(os.listdir(path))
    files_to_read = ROOT.std.vector('string')()
    for i in range(n_files):
        files_to_read.push_back(path + filenames[i])

    start_time = time.time()

    data = ROOT.RDataFrame("rootuple/CandidateTree", files_to_read)
    cut = data.Filter("candidate_charge == 0")\
          .Filter("candidate_cosAlpha > 0.99")\
          .Filter("candidate_vProb > 0.05")\
          .Filter("candidate_lxy / candidate_lxyErr > 3.0")\
          .Filter("ditrack_mass > 1.014")\
          .Filter("ditrack_mass < 1.024")\
          .Filter("candidate_vMass > 5.33")\
          .Filter("candidate_vMass < 5.40")
    np_array = cut.AsNumpy(["candidate_vMass"])

    return time.time() - start_time

def runtime_measure_mt(path, n_files, n_threads):
    ROOT.ROOT.DisableImplicitMT()

    if n_files == 0: return 0
    if n_threads == 0:
        return runtime_measure(path, n_files, False)

    ROOT.ROOT.EnableImplicitMT(n_threads)

    # Get paths to all the files to be read
    filenames = sorted(os.listdir(path))
    files_to_read = ROOT.std.vector('string')()
    for i in range(n_files):
        files_to_read.push_back(path + filenames[i])

    # Measure runtime
    start_time = time.time()

    data = ROOT.RDataFrame("rootuple/CandidateTree", files_to_read)
    cut = data.Filter("candidate_charge == 0")\
          .Filter("candidate_cosAlpha > 0.99")\
          .Filter("candidate_vProb > 0.05")\
          .Filter("candidate_lxy / candidate_lxyErr > 3.0")\
          .Filter("ditrack_mass > 1.014")\
          .Filter("ditrack_mass < 1.024")\
          .Filter("candidate_vMass > 5.33")\
          .Filter("candidate_vMass < 5.40")
    np_array = cut.AsNumpy(["candidate_vMass"])

    return time.time() - start_time

def runtime_measure_mp(path, n_files, n_processes):
    ROOT.ROOT.DisableImplicitMT()
    if n_files == 0: return 0
    if n_processes == 0: return runtime_measure(path, n_files, False)
    start_time = time.time()
    partitions = partition(path, n_files, n_processes)
    processes = []
    result = multiprocessing.Manager().list()
    for files in partitions:
        p = multiprocessing.Process(target=to_numpy, args=[files, result])
        p.start()
        processes.append(p)

    for p in processes:
        p.join()

    np.concatenate(tuple(result))
    runtime = time.time() - start_time

    return runtime

def runtime_vs_variable(path, target_dir, measure_function, variable, step, n_loops, var_max, constant=None):
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)

    result_path = ("%s/runtime_vs_%s_%d_%d_%d_%d.csv" % (target_dir, variable, constant, var_max, step, n_loops)) if constant else ("%s/runtime_vs_%s_%d_%d_%d.csv" % (target_dir, variable, var_max, step, n_loops))

    x = [a for a in range(0, var_max + step, step)]

    if not os.path.exists(result_path):
        with open(result_path, "w", newline="") as f:
            csv.writer(f).writerow(x)
    print(">Warming up. First measurement will be ignored.")
    y = [measure_function(*(path, i if "size" in variable else constant, constant if "size" in variable else i) if constant else (path, i, True)) for i in range(0, step, step)]
    for n in tqdm(range(n_loops)):
        with open(result_path, "r") as f:
            if sum(1 for row in csv.reader(f)) == n_loops + 1: break
        y = [measure_function(*(path, i if "size" in variable else constant, constant if "size" in variable else i) if constant else (path, i, True)) for i in range(0, var_max + step, step)]
        with open(result_path, "a", newline="") as f:
            csv.writer(f).writerow(y)

usefullvars = ["size","threads","size_mt"]

@click.command()
@click.option('--step', default=4, help='Scanning step.')
@click.option('--loops', default=20, help='No. of measurements to do.')
@click.option('--max', default=128, help='Variable max.')
@click.option('--name', prompt='Path', default = "/lustrehome/hdhoang2001/data/128_files/",
              help='File path.')
@click.option('--output', prompt='Output', default = "runtime_tests_rdf/128_files/",
              help='Output dir.')
@click.option('--variable', prompt='Variable to test', default = "size",
              help='Running variable (either size or threads).')
@click.option('--n_threads', prompt='Number of threads', default = 32,
              help='Number of threads to run. If variable==threads this will be useless.')
@click.option('--n_files', prompt='Number of files', default = 128,
              help='Number of files to load (used when testing multiple processes variable).')
def run(step,loops,max,name,output,variable,n_files,n_threads):

    if variable not in usefullvars:
        raise Exception("Sorry, only %s as variables. No %s here!"%(usefullvars,variable))

    if variable == "size":
        runtime_vs_variable(name, output, runtime_measure, variable, step, loops, max)
    elif variable == "threads":
        runtime_vs_variable(name, output, runtime_measure_mt, variable, step, loops, max, n_threads)
    elif variable == "size_mt":
        runtime_vs_variable(name, output, runtime_measure_mt, variable, step, loops, max, n_threads)


if __name__ == '__main__':
    run()
