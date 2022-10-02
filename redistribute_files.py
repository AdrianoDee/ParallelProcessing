import uproot
import os
import multiprocessing
import math
import pandas as pd

path = "../merged/"

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
    while not bool(slices) or slices[-1][-1] != (file_entries[-1]):
        slices.append([file_start, entry_start] + partition_helper(slice_entries, file_entries, file_start, entry_start))
        file_start = slices[-1][-2]
        entry_start = slices[-1][-1]
    return slices

def write_one_file(candidate_trees, candidate_slices, ups_trees, ups_slices, index, target_dir):
    candidate_data = []
    ups_data = []
    for i in range(candidate_slices[index][0], candidate_slices[index][2] + 1):
        candidate_data.append(candidate_trees[i].arrays(
            [key for key in candidate_trees[i].keys() if not key.endswith("_p4")],
            entry_start=candidate_slices[index][1] if i == candidate_slices[index][0] else None,
            entry_stop=candidate_slices[index][3] if i == candidate_slices[index][2] else None,
            library="pd"))
    for i in range(ups_slices[index][0], ups_slices[index][2] + 1):
        ups_data.append(ups_trees[i].arrays(
            [key for key in ups_trees[i].keys() if not key.endswith("_p4")],
            entry_start=ups_slices[index][1] if i == ups_slices[index][0] else None,
            entry_stop=ups_slices[index][3] if i == ups_slices[index][2] else None,
            library="pd"))
    file = uproot.recreate(target_dir + "/file" + str(index) + ".root")
    file.mkdir("rootuple")
    file["rootuple/CandidateTree"] = pd.concat(candidate_data)
    file["rootuple/UpsTree"] = pd.concat(ups_data)

def redistribute(path, n_files):
    target_dir = "../data/" + str(n_files) + "_files"
    os.mkdir(target_dir)
    candidate_trees = [uproot.open(path=path + filename + ":rootuple/CandidateTree", object_cache=None, array_cache=None) for filename in sorted(os.listdir(path))]
    candidate_slices = partition(candidate_trees, n_files)
    ups_trees = [uproot.open(path=path + filename + ":rootuple/UpsTree", object_cache=None, array_cache=None) for filename in sorted(os.listdir(path))]
    ups_slices = partition(ups_trees, n_files)
    result = multiprocessing.Manager().list()
    processes = []
    for i in range(n_files):
        p = multiprocessing.Process(target=write_one_file, args=[candidate_trees, candidate_slices, ups_trees, ups_slices, i, target_dir])
        p.start()
        processes.append(p)

    for p in processes:
        p.join()
    
# redistribute(path, 32)  
# redistribute(path, 64)
# redistribute(path, 128)
# redistribute(path, 1)
