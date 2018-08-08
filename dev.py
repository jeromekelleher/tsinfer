
import random
import os
import h5py
import zarr
import sys
import pandas as pd
import daiquiri
#import bsddb3
import time
import scipy
import pickle
import collections
import itertools
import tqdm
import shutil
import pprint
import numpy as np
import json

import matplotlib as mp
# Force matplotlib to not use any Xwindows backend.
mp.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

import tsinfer
import tsinfer.eval_util as eval_util
import msprime



def plot_breakpoints(ts, map_file, output_file):
    # Read in the recombination map using the read_hapmap engine,
    recomb_map = msprime.RecombinationMap.read_hapmap(map_file)

    # Now we get the positions and rates from the recombination
    # map and plot these using 500 bins.
    positions = np.array(recomb_map.get_positions()[1:])
    rates = np.array(recomb_map.get_rates()[1:])
    num_bins = 500
    v, bin_edges, _ = scipy.stats.binned_statistic(
        positions, rates, bins=num_bins)
    x = bin_edges[:-1][np.logical_not(np.isnan(v))]
    y = v[np.logical_not(np.isnan(v))]
    fig, ax1 = plt.subplots(figsize=(16, 6))
    ax1.plot(x, y, color="blue", label="Recombination rate")
    ax1.set_ylabel("Recombination rate")
    ax1.set_xlabel("Chromosome position")

    # Now plot the density of breakpoints along the chromosome
    breakpoints = np.array(list(ts.breakpoints()))
    ax2 = ax1.twinx()
    v, bin_edges = np.histogram(breakpoints, num_bins, density=True)
    ax2.plot(bin_edges[:-1], v, color="green", label="Breakpoint density")
    ax2.set_ylabel("Breakpoint density")
    ax2.set_xlim(1.5e7, 5.3e7)
    plt.legend()
    fig.savefig(output_file)


def make_errors(v, p):
    """
    For each sample an error occurs with probability p. Errors are generated by
    sampling values from the stationary distribution, that is, if we have an
    allele frequency of f, a 1 is emitted with probability f and a
    0 with probability 1 - f. Thus, there is a possibility that an 'error'
    will in fact result in the same value.
    """
    w = np.copy(v)
    if p > 0:
        m = v.shape[0]
        frequency = np.sum(v) / m
        # Randomly choose samples with probability p
        samples = np.where(np.random.random(m) < p)[0]
        # Generate observations from the stationary distribution.
        errors = (np.random.random(samples.shape[0]) < frequency).astype(int)
        w[samples] = errors
    return w


def generate_samples(ts, error_p):
    """
    Returns samples with a bits flipped with a specified probability.

    Rejects any variants that result in a fixed column.
    """
    S = np.zeros((ts.sample_size, ts.num_mutations), dtype=np.int8)
    for variant in ts.variants():
        done = False
        # Reject any columns that have no 1s or no zeros
        while not done:
            S[:, variant.index] = make_errors(variant.genotypes, error_p)
            s = np.sum(S[:, variant.index])
            done = 0 < s < ts.sample_size
    return S.T


def simple_augment_sites(sample_data, ts, **kwargs):
    edges = sorted(ts.edges(), key=lambda x: (x.child, x.left))
    last_edge = edges[0]
    srb_map = collections.defaultdict(list)
    for edge in edges[1:]:
        if edge.child == last_edge.child and edge.left == last_edge.right:
            if ts.node(edge.child).is_sample():
                key = edge.left, last_edge.parent, edge.parent
                srb_map[key].append(edge.child)
        last_edge = edge

    srbs = []
    for k, v in srb_map.items():
        if len(v) > 1:
            # print(k, "\t", v)
            srbs.append((k[0], v))
    srbs.sort()

    augmented_samples = tsinfer.SampleData(
        sequence_length=sample_data.sequence_length, **kwargs)
    position = sample_data.sites_position[:]
    srb_iter = iter(srbs)
    x, samples = next(srb_iter)
    for j, genotypes in sample_data.genotypes():
        # print(j, position[j], genotypes)
        inserted_genotypes = []
        while x == position[j]:
            a = np.zeros_like(genotypes)
            a[samples] = 1
            inserted_genotypes.append(a)
            x, samples = next(srb_iter, (-1, None))
        if len(inserted_genotypes) > 0:
            # Put the inserted genotypes before pos.
            distance = position[j] - position[j - 1]
            delta = distance / (len(inserted_genotypes) + 1)
            y = position[j - 1] + delta
            for inserted_genotype in inserted_genotypes:
                # print("\tAugmented site @ ", y)
                augmented_samples.add_site(position=y, genotypes=inserted_genotype)
                y += delta
        augmented_samples.add_site(position=position[j], genotypes=genotypes)
    augmented_samples.finalise()
    return augmented_samples

def augment_sites(sample_data, ts, **kwargs):
    # We're only interested in sample nodes for now.
    flags = ts.tables.nodes.flags
    is_sample = (flags & 1) != 0  # Only if bit 1 is set
    print(is_sample)

    edges = ts.tables.edges
    index = is_sample[edges.child]
    print(index)
    print(edges)

    left = edges.left
    right = edges.right
    parent = edges.parent
    child = edges.child

    order = np.lexsort([left, child])
    print(order)


    edges = sorted(ts.edges(), key=lambda x: (x.child, x.left))
    for j, edge in enumerate(edges):
        assert left[order[j]] == edge.left
        assert right[order[j]] == edge.right
        assert child[order[j]] == edge.child
        assert parent[order[j]] == edge.parent


    last_edge = edges[0]
    srb_map = collections.defaultdict(list)
    j = 1
    for edge in edges[1:]:
        if edge.child == last_edge.child and edge.left == last_edge.right:
            if ts.node(edge.child).is_sample():
                key = edge.left, last_edge.parent, edge.parent
                srb_map[key].append(edge.child)
        last_edge = edge
        j += 1



def tsinfer_dev(
        n, L, seed, num_threads=1, recombination_rate=1e-8,
        error_rate=0, engine="C", log_level="WARNING",
        debug=True, progress=False, path_compression=True):

    np.random.seed(seed)
    random.seed(seed)
    L_megabases = int(L * 10**6)

    # daiquiri.setup(level=log_level)

    source_ts = msprime.simulate(
            n, Ne=10**4, length=L_megabases,
            recombination_rate=recombination_rate, mutation_rate=1e-8,
            random_seed=seed)
    assert source_ts.num_sites > 0

    sample_data = tsinfer.SampleData.from_tree_sequence(source_ts)

    ancestor_data = tsinfer.generate_ancestors(
        sample_data, engine=engine, num_threads=num_threads)
    ancestors_ts = tsinfer.match_ancestors(
        sample_data, ancestor_data, engine=engine,
        path_compression=path_compression)

    ts = tsinfer.match_samples(sample_data, ancestors_ts,
            path_compression=path_compression, simplify=True, engine=engine)

    # augmented_samples = simple_augment_sites(sample_data, ts)
    augmented_samples = augment_sites(sample_data, ts)

#     edges = sorted(ts.edges(), key=lambda x: (x.child, x.left))
#     last_edge = edges[0]
#     srb_map = collections.defaultdict(list)
#     for edge in edges[1:]:
#         if edge.child == last_edge.child and edge.left == last_edge.right:
#             if ts.node(edge.child).is_sample():
#                 key = edge.left, last_edge.parent, edge.parent
#                 srb_map[key].append(edge.child)
#         last_edge = edge

#     srbs = []
#     for k, v in srb_map.items():
#         if len(v) > 1:
#             # print(k, "\t", v)
#             srbs.append((k[0], v))
#     srbs.sort()

#     augmented_samples = tsinfer.SampleData(sequence_length=sample_data.sequence_length)
#     position = sample_data.sites_position[:]
#     srb_iter = iter(srbs)
#     x, samples = next(srb_iter)
#     for j, genotypes in sample_data.genotypes():
#         # print(j, position[j], genotypes)
#         inserted_genotypes = []
#         while x == position[j]:
#             a = np.zeros_like(genotypes)
#             a[samples] = 1
#             inserted_genotypes.append(a)
#             x, samples = next(srb_iter, (-1, None))
#         if len(inserted_genotypes) > 0:
#             # Put the inserted genotypes before pos.
#             distance = position[j] - position[j - 1]
#             delta = distance / (len(inserted_genotypes) + 1)
#             y = position[j - 1] + delta
#             for inserted_genotype in inserted_genotypes:
#                 # print("\tAugmented site @ ", y)
#                 augmented_samples.add_site(position=y, genotypes=inserted_genotype)
#                 y += delta
#         augmented_samples.add_site(position=position[j], genotypes=genotypes)

#     augmented_samples.finalise()
    print("Sites", sample_data.num_sites)
    print("Added", augmented_samples.num_sites - sample_data.num_sites, "augmented sites")

    final_ts = tsinfer.infer(augmented_samples)
    print("nodes:", ts.num_nodes, final_ts.num_nodes, source_ts.num_nodes, sep="\t")
    print("edges:", ts.num_edges, final_ts.num_edges, source_ts.num_edges, sep="\t")
    print("trees:", ts.num_trees, final_ts.num_trees, source_ts.num_trees, sep="\t")
    sys.stdout.flush()

    breakpoints, kc_distance = eval_util.compare(ts, source_ts)
    d = breakpoints[1:] - breakpoints[:-1]
    d /= breakpoints[-1]
    no_augment = np.sum(kc_distance * d)
    breakpoints, kc_distance = eval_util.compare(final_ts, source_ts)
    d = breakpoints[1:] - breakpoints[:-1]
    d /= breakpoints[-1]
    augment = np.sum(kc_distance * d)
    print("kc   :", no_augment, augment)

#     # pos = ts.sequence_length / 8
#     pos = 0.001
#     for tree in ts.trees():
#         if tree.interval[0] <= pos < tree.interval[1]:
#             break
#     print(tree.draw(format="unicode"))

#     for tree in final_ts.trees():
#         if tree.interval[0] <= pos < tree.interval[1]:
#             break
#     print(tree.draw(format="unicode"))

#     for tree in source_ts.trees():
#         if tree.interval[0] <= pos < tree.interval[1]:
#             break
#     print(tree.draw(format="unicode"))


def dump_provenance(ts):
    print("dump provenance")
    for p in ts.provenances():
        print("-" * 50)
        print(p.timestamp)
        pprint.pprint(json.loads(p.record))


def build_profile_inputs(n, num_megabases):
    L = num_megabases * 10**6
    input_file = "tmp__NOBACKUP__/profile-n={}-m={}.input.trees".format(
            n, num_megabases)
    if os.path.exists(input_file):
        ts = msprime.load(input_file)
    else:
        ts = msprime.simulate(
            n, length=L, Ne=10**4, recombination_rate=1e-8, mutation_rate=1e-8,
            random_seed=10)
        print("Ran simulation: n = ", n, " num_sites = ", ts.num_sites,
                "num_trees =", ts.num_trees)
        ts.dump(input_file)
    filename = "tmp__NOBACKUP__/profile-n={}-m={}.samples".format(n, num_megabases)
    if os.path.exists(filename):
        os.unlink(filename)
    # daiquiri.setup(level="DEBUG")
    with tsinfer.SampleData(
            sequence_length=ts.sequence_length, path=filename,
            num_flush_threads=4) as sample_data:
        # progress_monitor = tqdm.tqdm(total=ts.num_samples)
        # for j in range(ts.num_samples):
        #     sample_data.add_sample(metadata={"name": "sample_{}".format(j)})
        #     progress_monitor.update()
        # progress_monitor.close()
        progress_monitor = tqdm.tqdm(total=ts.num_sites)
        for variant in ts.variants():
            sample_data.add_site(variant.site.position, variant.genotypes)
            progress_monitor.update()
        progress_monitor.close()

    print(sample_data)

#     filename = "tmp__NOBACKUP__/profile-n={}_m={}.ancestors".format(n, num_megabases)
#     if os.path.exists(filename):
#         os.unlink(filename)
#     ancestor_data = tsinfer.AncestorData.initialise(sample_data, filename=filename)
#     tsinfer.build_ancestors(sample_data, ancestor_data, progress=True)
#     ancestor_data.finalise()

def copy_1kg():
    source = "tmp__NOBACKUP__/1kg_chr22.samples"
    sample_data = tsinfer.SampleData.load(source)
    copy = sample_data.copy("tmp__NOBACKUP__/1kg_chr22_copy.samples")
    copy.finalise()
    print(sample_data)
    print("copy = ")
    print(copy)

def tutorial_samples():
    import tqdm
    import msprime
    import tsinfer

    ts = msprime.simulate(
        sample_size=10000, Ne=10**4, recombination_rate=1e-8,
        mutation_rate=1e-8, length=10*10**6, random_seed=42)
    ts.dump("tmp__NOBACKUP__/simulation-source.trees")
    print("simulation done:", ts.num_trees, "trees and", ts.num_sites,  "sites")

    progress = tqdm.tqdm(total=ts.num_sites)
    with tsinfer.SampleData(
            path="tmp__NOBACKUP__/simulation.samples",
            sequence_length=ts.sequence_length,
            num_flush_threads=2) as sample_data:
        for var in ts.variants():
            sample_data.add_site(var.site.position, var.genotypes, var.alleles)
            progress.update()
    progress.close()


def subset_sites(ts, position):
    """
    Return a copy of the specified tree sequence with sites reduced to those
    with positions in the specified list.
    """
    tables = ts.dump_tables()
    lookup = frozenset(position)
    tables.sites.clear()
    tables.mutations.clear()
    for site in ts.sites():
        if site.position in lookup:
            site_id = tables.sites.add_row(
                site.position, ancestral_state=site.ancestral_state,
                metadata=site.metadata)
            for mutation in site.mutations:
                tables.mutations.add_row(
                    site_id, node=mutation.node, parent=mutation.parent,
                    derived_state=mutation.derived_state,
                    metadata=mutation.metadata)
    return tables.tree_sequence()

def minimise(ts):
    tables = ts.dump_tables()

    out_map = {}
    in_map = {}
    first_site = 0
    for (_, edges_out, edges_in), tree in zip(ts.edge_diffs(), ts.trees()):
        for edge in edges_out:
            out_map[edge.child] = edge
        for edge in edges_in:
            in_map[edge.child] = edge
        if tree.num_sites > 0:
            sites = list(tree.sites())
            if first_site:
                x = 0
                first_site = False
            else:
                x = sites[0].position
            print("X = ", x)
            for edge in out_map.values():
                print("FLUSH", edge)
            for edge in in_map.values():
                print("INSER", edge)

            # # Flush the edge buffer.
            # for left, parent, child in edge_buffer:
            #     tables.edges.add_row(left, x, parent, child)
            # # Add edges for each node in the tree.
            # edge_buffer.clear()
            # for root in tree.roots:
            #     for u in tree.nodes(root):
            #         if u != root:
            #             edge_buffer.append((x, tree.parent(u), u))

    # position = np.hstack([[0], tables.sites.position, [ts.sequence_length]])
    # position = tables.sites.position
    # edges = []
    # print(position)
    # tables.edges.clear()
    # for edge in ts.edges():
    #     left = np.searchsorted(position, edge.left)
    #     right = np.searchsorted(position, edge.right)

    #     print(edge, left, right)
    #     # if right - left > 1:
    #         # print("KEEP:", edge, left, right)
    #         # tables.edges.add_row(
    #         #     position[left], position[right], edge.parent, edge.child)
    #         # print("added", tables.edges[-1])
    #     # else:
    #         # print("SKIP:", edge, left, right)

    # ts = tables.tree_sequence()
    # for tree in ts.trees():
    #     print("TREE:", tree.interval)
    #     print(tree.draw(format="unicode"))


def run_build():

    sample_data = tsinfer.load(sys.argv[1])
    ad = tsinfer.generate_ancestors(sample_data)
    print(ad)


if __name__ == "__main__":

    # run_build()

    # np.set_printoptions(linewidth=20000)
    # np.set_printoptions(threshold=20000000)

    # tutorial_samples()

    # build_profile_inputs(10, 10)
    # build_profile_inputs(100, 10)
    # build_profile_inputs(1000, 100)
    # build_profile_inputs(10**4, 100)
    # build_profile_inputs(10**5, 100)

    # for j in range(1, 100):
    #     tsinfer_dev(15, 0.5, seed=j, num_threads=0, engine="P", recombination_rate=1e-8)
    # copy_1kg()
    # tsinfer_dev(105, 10.25, seed=4, num_threads=0, engine="C", recombination_rate=2e-8,
    #         path_compression=True)
    tsinfer_dev(5, 0.1, seed=4, num_threads=0, engine="C", recombination_rate=1e-8,
            path_compression=True)


#     for seed in range(1, 10000):
#         print(seed)
#         # tsinfer_dev(40, 2.5, seed=seed, num_threads=1, genotype_quality=1e-3, engine="C")
