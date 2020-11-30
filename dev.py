import random
import os
import sys

import tqdm
import pprint
import numpy as np
import json
import msprime

# import tskit

import tsinfer


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


def tsinfer_dev(
    n,
    L,
    seed,
    num_threads=1,
    recombination_rate=1e-8,
    error_rate=0,
    engine="C",
    log_level="WARNING",
    precision=None,
    debug=True,
    progress=False,
    path_compression=True,
):

    np.random.seed(seed)
    random.seed(seed)
    L_megabases = int(L * 10 ** 6)

    # daiquiri.setup(level=log_level)

    source_ts = msprime.simulate(
        n,
        Ne=10 ** 4,
        length=L_megabases,
        recombination_rate=recombination_rate,
        mutation_rate=1e-8,
        random_seed=seed,
    )
    if debug:
        print("num_sites = ", source_ts.num_sites)
    assert source_ts.num_sites > 0

    # ts = msprime.mutate(ts, rate=1e-8, random_seed=seed,
    #         model=msprime.InfiniteSites(msprime.NUCLEOTIDES))

    # samples = tsinfer.SampleData.from_tree_sequence(ts)

    with tsinfer.SampleData(sequence_length=source_ts.sequence_length) as samples:
        for var in source_ts.variants():
            # var.genotypes[var.site.id % source_ts.num_samples] = tskit.MISSING_DATA
            samples.add_site(var.site.position, var.genotypes, var.alleles)

    print(samples)
    # for variant in samples.variants():
    #     print(variant)

    rho = recombination_rate
    mu = 1e-3  # 1e-15

    #     num_alleles = samples.num_alleles(inference_sites=True)
    #     num_sites = samples.num_inference_sites
    #     with tsinfer.AncestorData(samples) as ancestor_data:
    #         t = np.sum(num_alleles) + 1
    #         for j in range(num_sites):
    #             for allele in range(num_alleles[j]):
    #                 ancestor_data.add_ancestor(j, j + 1, t, [j], [allele])
    #                 t -= 1

    ancestor_data = tsinfer.generate_ancestors(
        samples, engine=engine, num_threads=num_threads
    )
    print(ancestor_data)

    ancestors_ts = tsinfer.match_ancestors(
        samples,
        ancestor_data,
        engine=engine,
        path_compression=True,
        extended_checks=False,
        precision=precision,
        recombination_rate=rho,
        mismatch_rate=mu,
    )
    # print(ancestors_ts.tables)

    # print("ancestors ts")
    # for tree in ancestors_ts.trees():
    #     print(tree.draw_text())
    #     for site in tree.sites():
    #         if len(site.mutations) > 1:
    #             print(site.id)
    #             for mutation in site.mutations:
    #                 print("\t", mutation.node, mutation.derived_state)

    # for var in ancestors_ts.variants():
    #     print(var.genotypes)

    # print(ancestors_ts.tables)

    # ancestors_ts = tsinfer.augment_ancestors(samples, ancestors_ts,
    #         [5, 6, 7], engine=engine)

    ts = tsinfer.match_samples(
        samples,
        ancestors_ts,
        recombination_rate=rho,
        mismatch_rate=mu,
        path_compression=False,
        engine=engine,
        precision=precision,
        simplify=False,
    )
    # print(ts.tables)

    #     for var1, var2, var3 in zip(
    #         source_ts.variants(), ts.variants(), samples.variants()
    #     ):
    #         if np.any(var1.genotypes != var2.genotypes):
    #             print("mismatch at ", var1.site.id)
    #             print(var1.genotypes)
    #             print(var2.genotypes)
    #             print(var3.genotypes)

    print("num_edges = ", ts.num_edges)

    # # print(ts.draw_text())
    # for tree in ts.trees():
    #     print(tree.draw_text())
    #     for site in tree.sites():
    #         if len(site.mutations) > 1:
    #             print(site.id)
    #             for mutation in site.mutations:
    #                 print("\t", mutation.node, mutation.derived_state)

    # # print(ts.tables.edges)
    # print(ts.dump_tables())

    # simplified = ts.simplify()
    # print("edges before = ", simplified.num_edges)

    # new_ancestors_ts = insert_srb_ancestors(ts)
    # ts = tsinfer.match_samples(samples, new_ancestors_ts,
    #         path_compression=False, engine=engine,
    #         simplify=True)

    #     for tree in ts.trees():
    #         print(tree.interval)
    #         print(tree.draw(format="unicode"))

    # print(ts.tables.edges)
    # for tree in ts.trees():
    #     print(tree.draw(format="unicode"))

    tsinfer.verify(samples, ts)


#     # output_ts = tsinfer.match_samples(subset_samples, ancestors_ts, engine=engine)
#     output_ts = tsinfer.match_samples(sample_data, ancestors_ts, engine=engine)
#     # dump_provenance(output_ts)


def dump_provenance(ts):
    print("dump provenance")
    for p in ts.provenances():
        print("-" * 50)
        print(p.timestamp)
        pprint.pprint(json.loads(p.record))


def build_profile_inputs(n, num_megabases):
    L = num_megabases * 10 ** 6
    input_file = "tmp__NOBACKUP__/profile-n={}-m={}.input.trees".format(
        n, num_megabases
    )
    if os.path.exists(input_file):
        ts = msprime.load(input_file)
    else:
        ts = msprime.simulate(
            n,
            length=L,
            Ne=10 ** 4,
            recombination_rate=1e-8,
            mutation_rate=1e-8,
            random_seed=10,
        )
        print(
            "Ran simulation: n = ",
            n,
            " num_sites = ",
            ts.num_sites,
            "num_trees =",
            ts.num_trees,
        )
        ts.dump(input_file)
    filename = "tmp__NOBACKUP__/profile-n={}-m={}.samples".format(n, num_megabases)
    if os.path.exists(filename):
        os.unlink(filename)
    # daiquiri.setup(level="DEBUG")
    with tsinfer.SampleData(
        sequence_length=ts.sequence_length, path=filename, num_flush_threads=4
    ) as sample_data:
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
        sample_size=10000,
        Ne=10 ** 4,
        recombination_rate=1e-8,
        mutation_rate=1e-8,
        length=10 * 10 ** 6,
        random_seed=42,
    )
    ts.dump("tmp__NOBACKUP__/simulation-source.trees")
    print("simulation done:", ts.num_trees, "trees and", ts.num_sites, "sites")

    progress = tqdm.tqdm(total=ts.num_sites)
    with tsinfer.SampleData(
        path="tmp__NOBACKUP__/simulation.samples",
        sequence_length=ts.sequence_length,
        num_flush_threads=2,
    ) as sample_data:
        for var in ts.variants():
            sample_data.add_site(var.site.position, var.genotypes, var.alleles)
            progress.update()
    progress.close()


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

    # for j in range(1, 1000):
    # for j in [96]:
    #     print("seed = ", j)
    #     tsinfer_dev(5, 0.5, seed=j, num_threads=0, engine="P", recombination_rate=1e-8)
    # copy_1kg()
    tsinfer_dev(
        8,
        0.05,
        seed=4,
        num_threads=0,
        engine="C",
        recombination_rate=1e-8,
        precision=0,
    )

    # minimise_dev()

#     for seed in range(1, 10000):
#         print(seed)
#         tsinfer_dev(40, 2.5, seed=seed, num_threads=1, genotype_quality=1e-3,
#         engine="C")
