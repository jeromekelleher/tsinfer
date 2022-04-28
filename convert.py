import datetime
import sys

import cyvcf2
import numpy as np
import pandas as pd
import tqdm
import tskit

import tsinfer


# https://stackoverflow.com/questions/6451655/
def year_fraction(datestr):
    # if len(datestr.split("-")) == 2:
    # datestr += "-01"
    date = datetime.date.fromisoformat(datestr)
    start = datetime.date(date.year, 1, 1).toordinal()
    year_length = datetime.date(date.year + 1, 1, 1).toordinal() - start
    return date.year + float(date.toordinal() - start) / year_length


def add_sites(vcf, sample_data, keep_samples):
    pbar = tqdm.tqdm(total=sample_data.sequence_length)
    pos = 0
    for variant in vcf:  # Loop over variants, each assumed at a unique site
        pbar.update(variant.POS - pos)
        if pos == variant.POS:
            raise ValueError("Duplicate positions for variant at position", pos)
        else:
            pos = variant.POS
        # print(pos, samples.sequence_length)
        # Assume REF is the ancestral state.
        alleles = [variant.REF] + variant.ALT
        genotypes = np.array(variant.genotypes).T[0]
        sample_data.add_site(pos, genotypes=genotypes[keep_samples], alleles=alleles)
    pbar.close()


def convert_vcf():

    vcf = cyvcf2.VCF("merged.vcf.gz")
    df = pd.read_csv("merged.tsv", sep="\t")
    # print(df)
    df.Location.fillna("NA", inplace=True)
    df.Collection_Date.fillna("NA", inplace=True)
    sample_metadata_map = {}
    for _, row in df.iterrows():
        # print(row)
        values = {
            "virus": row["Virus"],
            "collection_date": row["Collection_Date"],
            "location": row["Location"],
        }
        sample_metadata_map[row["Accession_ID"]] = values
        for iid in row["Identical_Seq"]:
            sample_metadata_map[iid] = values

    print("loaded metadata")

    with tsinfer.SampleData(path="merged.samples", sequence_length=29904) as samples:
        population_map = {}
        for population in sorted(set(df["Location"])):
            if population != "NA":
                pid = samples.add_population(metadata={"name": population})
            else:
                pid = -1
            population_map[population] = pid

        keep_samples = np.ones(len(vcf.samples), dtype=bool)
        for j, sample in enumerate(vcf.samples):
            metadata = None
            try:
                metadata = dict(sample_metadata_map[sample])
                time = year_fraction(metadata["collection_date"])
                population = population_map[metadata.pop("location")]
                metadata["name"] = sample
                samples.add_individual(
                    ploidy=1, population=population, time=time, metadata=metadata
                )
            except (KeyError, ValueError):
                # print("Discarding", sample, metadata)
                keep_samples[j] = 0

        print("Discarded", np.sum(keep_samples == 0))

        add_sites(vcf, samples, keep_samples)


def naive_inference():
    sd = tsinfer.load("merged-trimmed.samples")
    time = sd.individuals_time[:]

    present_zero = np.max(time) - time

    with sd.copy() as sd_copy:
        sd_copy.individuals_time[:] = present_zero

    a = np.zeros(sd.num_sites, dtype=int)
    # FIXME: this second ancestor shouldn't be necessary.
    with tsinfer.AncestorData(sd_copy) as ad:
        ad.add_ancestor(start=0, end=ad.num_sites, time=2, focal_sites=[], haplotype=a)
        ad.add_ancestor(start=0, end=ad.num_sites, time=1, focal_sites=[], haplotype=a)

    ancestors_ts = tsinfer.match_ancestors(sd_copy, ad)
    # print(ancestors_ts)

    pm = tsinfer.inference._get_progress_monitor(
        True, generate_ancestors=False, match_ancestors=False, match_samples=True,
    )

    for t in np.unique(present_zero)[::-1]:
        ids = list(np.where(present_zero == t)[0])
        print("t = ", t, ":", len(ids))
        sdt = sd_copy.subset(individuals=ids)
        # print(sdt)
        ts = tsinfer.match_samples(
            sdt,
            ancestors_ts,
            recombination_rate=1e-20,
            mismatch_ratio=1e10,
            num_threads=40,
            simplify=False,
            progress_monitor=True
        )
        print(ts)
        # print(ts.draw_text())
        # print(ts.tables.nodes)
        ts.dump(f"tmp__NOBACKUP__/covid/norecomb/t={t}.trees")
        tables = ts.dump_tables()
        tables.nodes.time += 1
        tables.nodes.population = np.zeros_like(tables.nodes.population) - 1
        tables.populations.clear()
        ancestors_ts = tables.tree_sequence()


def fixup_samples_metadata(ts):
    print(ts)
    tables = ts.dump_tables()
    nt = tables.nodes
    nt.truncate(2)
    nodes = iter(ts.nodes())
    for _ in range(2):
        next(nodes)
    for node in nodes:
        ind = ts.individual(node.individual)
        nt.add_row(flags=1, time=node.time, metadata=ind.metadata)
        # break
    nt.metadata_schema = tskit.MetadataSchema.permissive_json()
    tables.populations.clear()
    tables.individuals.clear()
    return tables.tree_sequence()

    # for node in ts.nodes():

    # print(tables.nodes)
    # tables.popaltions.clear()

def trim_data():
    source = tsinfer.load("merged.samples")
    print(source)

    keep_sites = []
    for variant in source.variants():
        if 100 < variant.site.position < 29000:
            is_snp = all(len(allele) < 5 for allele in variant.alleles)
            if is_snp:
                keep_sites.append(variant.site.id)

    print("running subset", len(keep_sites))
    source.subset(sites=keep_sites, path="merged-trimmed.samples", num_flush_threads=8)



if __name__ == "__main__":

    # print(df)
    # convert_vcf()
    # naive_inference()

    ts = fixup_samples_metadata(tskit.load(sys.argv[1]))
    ts.dump(sys.argv[2])
    # trim_data()
