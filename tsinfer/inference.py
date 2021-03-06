#
# Copyright (C) 2018 University of Oxford
#
# This file is part of tsinfer.
#
# tsinfer is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# tsinfer is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with tsinfer.  If not, see <http://www.gnu.org/licenses/>.
#
"""
TODO module docs.
"""

import collections
import queue
import time
import pickle
import logging
import threading

import numpy as np
import tqdm
import humanize
import msprime

import _tsinfer
import tsinfer.formats as formats
import tsinfer.algorithm as algorithm
import tsinfer.threads as threads

logger = logging.getLogger(__name__)

UNKNOWN_ALLELE = 255


# TODO figure out what this function is really for, and how we should
# parameterise it. At the moment it's just used as a convenience for
# round trip testing.
def infer(
        genotypes, positions=None, sequence_length=None,
        method="C", num_threads=0, progress=False, path_compression=True):

    num_sites, num_samples = genotypes.shape
    sample_data = formats.SampleData.initialise(
        num_samples=num_samples, sequence_length=sequence_length, compressor=None)
    for j in range(num_sites):
        pos = j
        if positions is not None:
            pos = positions[j]
        sample_data.add_variant(pos, ["0", "1"], genotypes[j])
    sample_data.finalise()

    ancestor_data = formats.AncestorData.initialise(sample_data, compressor=None)
    build_ancestors(sample_data, ancestor_data, method=method)
    ancestor_data.finalise()

    ancestors_ts = match_ancestors(
        sample_data, ancestor_data, method=method, num_threads=num_threads,
        path_compression=path_compression)

    inferred_ts = match_samples(
        sample_data, ancestors_ts, method=method, num_threads=num_threads,
        path_compression=path_compression)
    return inferred_ts


def build_ancestors(input_data, ancestor_data, progress=False, method="C"):

    num_sites = input_data.num_variant_sites
    num_samples = input_data.num_samples
    if method == "C":
        logger.debug("Using C AncestorBuilder implementation")
        ancestor_builder = _tsinfer.AncestorBuilder(num_samples, num_sites)
    else:
        logger.debug("Using Python AncestorBuilder implementation")
        ancestor_builder = algorithm.AncestorBuilder(num_samples, num_sites)

    progress_monitor = tqdm.tqdm(total=num_sites, disable=not progress)
    frequency = input_data.frequency[:]
    logger.info("Starting site addition")
    for j, (site_id, genotypes) in enumerate(input_data.variants()):
        ancestor_builder.add_site(j, int(frequency[site_id]), genotypes)
        progress_monitor.update()
    progress_monitor.close()
    logger.info("Finished adding sites")

    descriptors = ancestor_builder.ancestor_descriptors()
    if len(descriptors) > 0:
        num_ancestors = len(descriptors)
        logger.info("Starting build for {} ancestors".format(num_ancestors))
        a = np.zeros(num_sites, dtype=np.uint8)
        root_time = descriptors[0][0] + 1
        ultimate_ancestor_time = root_time + 1
        # Add the ultimate ancestor. This is an awkward hack really; we don't
        # ever insert this ancestor. The only reason to add it here is that
        # it makes sure that the ancestor IDs we have in the ancestor file are
        # the same as in the ancestor tree sequence. This seems worthwhile.
        ancestor_data.add_ancestor(
            start=0, end=num_sites, time=ultimate_ancestor_time,
            focal_sites=[], haplotype=a)
        # Hack to ensure we always have a root with zeros at every position.
        ancestor_data.add_ancestor(
            start=0, end=num_sites, time=root_time,
            focal_sites=np.array([], dtype=np.int32), haplotype=a)
        progress_monitor = tqdm.tqdm(total=len(descriptors), disable=not progress)
        for freq, focal_sites in descriptors:
            before = time.perf_counter()
            # TODO: This is a read-only process so we can multithread it.
            s, e = ancestor_builder.make_ancestor(focal_sites, a)
            assert np.all(a[s: e] != UNKNOWN_ALLELE)
            assert np.all(a[:s] == UNKNOWN_ALLELE)
            assert np.all(a[e:] == UNKNOWN_ALLELE)
            duration = time.perf_counter() - before
            logger.debug(
                "Made ancestor with {} focal sites and length={} in {:.2f}s.".format(
                    focal_sites.shape[0], e - s, duration))
            ancestor_data.add_ancestor(
                start=s, end=e, time=freq, focal_sites=focal_sites, haplotype=a)
            progress_monitor.update()
        progress_monitor.close()
    logger.info("Finished building ancestors")


def match_ancestors(
        sample_data, ancestor_data, output_path=None, method="C", progress=False,
        num_threads=0, path_compression=True, traceback_file_pattern=None,
        extended_checks=False):
    """
    Runs the copying process of the specified input and ancestors and returns
    the resulting tree sequence.
    """
    matcher = AncestorMatcher(
        sample_data, ancestor_data, output_path=output_path, method=method,
        progress=progress, path_compression=path_compression,
        num_threads=num_threads, traceback_file_pattern=traceback_file_pattern,
        extended_checks=extended_checks)
    return matcher.match_ancestors()


def verify(input_hdf5, ancestors_hdf5, ancestors_ts, progress=False):
    """
    Runs the copying process of the specified input and ancestors and returns
    the resulting tree sequence.
    """
    input_file = formats.InputFile(input_hdf5)
    ancestor_data = formats.AncestorFile(ancestors_hdf5, input_file, 'r')
    # TODO change these value errors to VerificationErrors or something.
    if ancestors_ts.num_nodes != ancestor_data.num_ancestors:
        raise ValueError("Incorrect number of ancestors")
    if ancestors_ts.num_sites != input_file.num_sites:
        raise ValueError("Incorrect number of sites")

    progress_monitor = tqdm.tqdm(
        total=ancestors_ts.num_sites, disable=not progress, dynamic_ncols=True)

    count = 0
    for g1, v in zip(ancestor_data.site_genotypes(), ancestors_ts.variants()):
        g2 = v.genotypes
        # Set anything unknown to 0
        g1[g1 == UNKNOWN_ALLELE] = 0
        if not np.array_equal(g1, g2):
            raise ValueError("Unequal genotypes at site", v.id)
        progress_monitor.update()
        count += 1
    if count != ancestors_ts.num_sites:
        raise ValueError("Iteration stopped early")
    progress_monitor.close()


def match_samples(
        sample_data, ancestors_ts, method="C", progress=False,
        num_threads=0, path_compression=True, simplify=True,
        traceback_file_pattern=None, extended_checks=False,
        stabilise_node_ordering=False):
    manager = SampleMatcher(
        sample_data, ancestors_ts, path_compression=path_compression,
        method=method, progress=progress, num_threads=num_threads,
        traceback_file_pattern=traceback_file_pattern,
        extended_checks=extended_checks)
    manager.match_samples()
    ts = manager.finalise(
        simplify=simplify, stabilise_node_ordering=stabilise_node_ordering)
    return ts


class Matcher(object):

    # The description for the progress monitor bar.
    progress_bar_description = None

    def __init__(
            self, sample_data, num_threads=1, method="C",
            path_compression=True, progress=False, traceback_file_pattern=None,
            extended_checks=False):
        self.sample_data = sample_data
        self.num_threads = num_threads
        self.path_compression = path_compression
        self.num_samples = self.sample_data.num_samples
        self.num_sites = self.sample_data.num_variant_sites
        self.progress = progress
        self.extended_checks = extended_checks
        self.tree_sequence_builder_class = algorithm.TreeSequenceBuilder

        if method == "C":
            logger.debug("Using C matcher implementation")
            self.tree_sequence_builder_class = _tsinfer.TreeSequenceBuilder
            self.ancestor_matcher_class = _tsinfer.AncestorMatcher
        elif method == "Py-matrix":
            logger.debug("Using Python matrix implementation")
            self.ancestor_matcher_class = algorithm.MatrixAncestorMatcher
        else:
            logger.debug("Using Python matcher implementation")
            self.ancestor_matcher_class = algorithm.AncestorMatcher
        self.tree_sequence_builder = None
        # Debugging. Set this to a file path like "traceback_{}.pkl" to store the
        # the tracebacks for each node ID and other debugging information.
        self.traceback_file_pattern = traceback_file_pattern

        # Allocate 64K nodes and edges initially. This will double as needed and will
        # quickly be big enough even for very large instances.
        max_edges = 64 * 1024
        max_nodes = 64 * 1024
        self.tree_sequence_builder = self.tree_sequence_builder_class(
            num_sites=self.num_sites, max_nodes=max_nodes, max_edges=max_edges)
        logger.debug("Allocated tree sequence builder with max_nodes={}".format(
            max_nodes))

        # Allocate the matchers and statistics arrays.
        num_threads = max(1, self.num_threads)
        self.match = [np.zeros(self.num_sites, np.uint8) for _ in range(num_threads)]
        self.results = ResultBuffer()
        self.mean_traceback_size = np.zeros(num_threads)
        self.num_matches = np.zeros(num_threads)
        self.matcher = [
            self.ancestor_matcher_class(
                self.tree_sequence_builder, extended_checks=self.extended_checks)
            for _ in range(num_threads)]
        # The progress monitor is allocated later by subclasses.
        self.progress_monitor = None

    def allocate_progress_monitor(self, total, initial=0, postfix=None):
        bar_format = (
            "{desc}{percentage:3.0f}%|{bar}"
            "| {n_fmt}/{total_fmt} [{elapsed}, {rate_fmt}{postfix}]")
        self.progress_monitor = tqdm.tqdm(
            desc=self.progress_bar_description, bar_format=bar_format,
            total=total, disable=not self.progress, initial=initial,
            smoothing=0.01, postfix=postfix, dynamic_ncols=True)

    def _find_path(self, child_id, haplotype, start, end, thread_index=0):
        """
        Finds the path of the specified haplotype and upates the results
        for the specified thread_index.
        """
        matcher = self.matcher[thread_index]
        match = self.match[thread_index]
        # print("Find path", child_id)
        left, right, parent = matcher.find_path(haplotype, start, end, match)
        # print("Done")
        self.results.set_path(child_id, left, right, parent)
        self.progress_monitor.update()
        self.mean_traceback_size[thread_index] += matcher.mean_traceback_size
        self.num_matches[thread_index] += 1
        logger.debug("matched node {}; num_edges={} tb_size={:.2f} match_mem={}".format(
            child_id, left.shape[0], matcher.mean_traceback_size,
            humanize.naturalsize(matcher.total_memory, binary=True)))
        if self.traceback_file_pattern is not None:
            # Write out the traceback debug. WARNING: this will be huge!
            filename = self.traceback_file_pattern.format(child_id)
            traceback = [matcher.get_traceback(l) for l in range(self.num_sites)]
            with open(filename, "wb") as f:
                debug = {
                    "child_id:": child_id,
                    "haplotype": haplotype,
                    "start": start,
                    "end": end,
                    "match": match,
                    "traceback": traceback}
                pickle.dump(debug, f)
                logger.debug(
                    "Dumped ancestor traceback debug to {}".format(filename))
        return left, right, parent

    def restore_tree_sequence_builder(self, ancestors_ts):
        # before = time.perf_counter()
        tables = ancestors_ts.dump_tables()
        nodes = tables.nodes
        self.tree_sequence_builder.restore_nodes(nodes.time, nodes.flags)
        edges = tables.edges
        # Need to sort by child ID here and left so that we can efficiently
        # insert the child paths.
        # TODO remove this step when we use a native zarr file for storing the
        # ancestor tree sequence. We output the edges in this order and we're
        # just sorting/resorting the edges here.
        index = np.lexsort((edges.left, edges.child))
        self.tree_sequence_builder.restore_edges(
            edges.left.astype(np.int32)[index],
            edges.right.astype(np.int32)[index],
            edges.parent[index],
            edges.child[index])
        mutations = tables.mutations
        self.tree_sequence_builder.restore_mutations(
            mutations.site, mutations.node, mutations.derived_state - ord('0'),
            mutations.parent)
        self.mutated_sites = mutations.site
        # print("SITE  =", self.mutated_sites)
        logger.info(
            "Loaded {} samples {} nodes; {} edges; {} sites; {} mutations".format(
                ancestors_ts.num_samples, len(nodes), len(edges), ancestors_ts.num_sites,
                len(mutations)))

    def get_tree_sequence(self, rescale_positions=True, all_sites=False):
        """
        Returns the current state of the build tree sequence. All samples and
        ancestors will have the sample node flag set.
        """
        # TODO Change the API here to ask whether we want a final tree sequence
        # or not. In the latter case we also need to translate the ancestral
        # and derived states to the input values.
        tsb = self.tree_sequence_builder
        flags, time = tsb.dump_nodes()
        nodes = msprime.NodeTable()
        nodes.set_columns(flags=flags, time=time)

        left, right, parent, child = tsb.dump_edges()
        if rescale_positions:
            position = self.sample_data.position[:]
            sequence_length = self.sample_data.sequence_length
            if sequence_length is None or sequence_length < position[-1]:
                sequence_length = position[-1] + 1
            # Subset down to the variants.
            position = position[self.sample_data.variant_site[:]]
            x = np.hstack([position, [sequence_length]])
            x[0] = 0
            left = x[left]
            right = x[right]
        else:
            position = np.arange(tsb.num_sites)
            sequence_length = max(1, tsb.num_sites)

        edges = msprime.EdgeTable()
        edges.set_columns(left=left, right=right, parent=parent, child=child)

        sites = msprime.SiteTable()
        sites.set_columns(
            position=position,
            ancestral_state=np.zeros(tsb.num_sites, dtype=np.int8) + ord('0'),
            ancestral_state_offset=np.arange(tsb.num_sites + 1, dtype=np.uint32))
        mutations = msprime.MutationTable()
        site = np.zeros(tsb.num_mutations, dtype=np.int32)
        node = np.zeros(tsb.num_mutations, dtype=np.int32)
        parent = np.zeros(tsb.num_mutations, dtype=np.int32)
        derived_state = np.zeros(tsb.num_mutations, dtype=np.int8)
        site, node, derived_state, parent = tsb.dump_mutations()
        derived_state += ord('0')
        mutations.set_columns(
            site=site, node=node, derived_state=derived_state,
            derived_state_offset=np.arange(tsb.num_mutations + 1, dtype=np.uint32),
            parent=parent)
        if all_sites:
            # Append the sites and mutations for each singleton.
            num_singletons = self.sample_data.num_singleton_sites
            singleton_site = self.sample_data.singleton_site[:]
            singleton_sample = self.sample_data.singleton_sample[:]
            pos = self.sample_data.position[:]
            new_sites = np.arange(
                len(sites), len(sites) + num_singletons, dtype=np.int32)
            sites.append_columns(
                position=pos[singleton_site],
                ancestral_state=np.zeros(num_singletons, dtype=np.int8) + ord('0'),
                ancestral_state_offset=np.arange(num_singletons + 1, dtype=np.uint32))
            mutations.append_columns(
                site=new_sites, node=self.sample_ids[singleton_sample],
                derived_state=np.zeros(num_singletons, dtype=np.int8) + ord('1'),
                derived_state_offset=np.arange(num_singletons + 1, dtype=np.uint32))
            # Get the invariant sites
            num_invariants = self.sample_data.num_invariant_sites
            invariant_site = self.sample_data.invariant_site[:]
            sites.append_columns(
                position=pos[invariant_site],
                ancestral_state=np.zeros(num_invariants, dtype=np.int8) + ord('0'),
                ancestral_state_offset=np.arange(num_invariants + 1, dtype=np.uint32))

        msprime.sort_tables(nodes, edges, sites=sites, mutations=mutations)
        return msprime.load_tables(
            nodes=nodes, edges=edges, sites=sites, mutations=mutations,
            sequence_length=sequence_length)


class AncestorMatcher(Matcher):
    progress_bar_description = "match-ancestors"

    def __init__(self, sample_data, ancestor_data, output_path, **kwargs):
        super().__init__(sample_data, **kwargs)
        self.output_path = output_path
        self.last_output_time = time.time()
        self.ancestor_data = ancestor_data
        self.num_ancestors = self.ancestor_data.num_ancestors
        self.epoch = self.ancestor_data.time[:]
        self.focal_sites = self.ancestor_data.focal_sites[:]
        self.start = self.ancestor_data.start[:]
        self.end = self.ancestor_data.end[:]

        # Create a list of all ID ranges in each epoch.
        if self.start.shape[0] == 0:
            self.num_epochs = 0
        else:
            breaks = np.where(self.epoch[1:] != self.epoch[:-1])[0]
            start = np.hstack([[0], breaks + 1])
            end = np.hstack([breaks + 1, [self.num_ancestors]])
            self.epoch_slices = np.vstack([start, end]).T
            self.num_epochs = self.epoch_slices.shape[0]
        first_ancestor = 1
        self.start_epoch = 1
        # Add nodes for all the ancestors so that the ancestor IDs are equal
        # to the node IDs.
        for ancestor_id in range(self.num_ancestors):
            self.tree_sequence_builder.add_node(self.epoch[ancestor_id])

        self.ancestors = self.ancestor_data.ancestors()
        if self.num_epochs > 0:
            # Consume the first ancestor.
            a = next(self.ancestors)
            assert np.array_equal(a, np.zeros(self.num_sites, dtype=np.uint8))
            self.allocate_progress_monitor(
                self.num_ancestors, initial=first_ancestor,
                postfix=self.__epoch_info_dict(self.start_epoch - 1))

    def __epoch_info_dict(self, epoch_index):
        start, end = self.epoch_slices[epoch_index]
        return collections.OrderedDict([
            ("edges", "{:.0G}".format(self.tree_sequence_builder.num_edges)),
            ("epoch", str(self.epoch[start])),
            ("nanc", str(end - start))
        ])

    def __update_progress_epoch(self, epoch_index):
        """
        Updates the progress monitor to show information about the present epoch
        """
        self.progress_monitor.set_postfix(self.__epoch_info_dict(epoch_index))

    def __ancestor_find_path(self, ancestor_id, ancestor, thread_index=0):
        haplotype = np.zeros(self.num_sites, dtype=np.uint8) + UNKNOWN_ALLELE
        focal_sites = self.focal_sites[ancestor_id]
        start = self.start[ancestor_id]
        end = self.end[ancestor_id]
        self.results.set_mutations(ancestor_id, focal_sites)
        assert ancestor.shape[0] == (end - start)
        haplotype[start: end] = ancestor
        assert np.all(haplotype[0: start] == UNKNOWN_ALLELE)
        assert np.all(haplotype[end:] == UNKNOWN_ALLELE)
        assert np.all(haplotype[focal_sites] == 1)
        logger.debug(
            "Finding path for ancestor {}; start={} end={} num_focal_sites={}".format(
                ancestor_id, start, end, focal_sites.shape[0]))
        haplotype[focal_sites] = 0
        left, right, parent = self._find_path(
                ancestor_id, haplotype, start, end, thread_index)
        assert np.all(self.match[thread_index] == haplotype)

    def __complete_epoch(self, epoch_index):
        start, end = map(int, self.epoch_slices[epoch_index])
        num_ancestors_in_epoch = end - start
        current_time = self.epoch[start]
        nodes_before = self.tree_sequence_builder.num_nodes

        for child_id in range(start, end):
            # TODO we should be adding the ancestor ID here as well as metadata.
            left, right, parent = self.results.get_path(child_id)
            self.tree_sequence_builder.add_path(
                child_id, left, right, parent, compress=self.path_compression)
            site, derived_state = self.results.get_mutations(child_id)
            self.tree_sequence_builder.add_mutations(child_id, site, derived_state)

        extra_nodes = (
            self.tree_sequence_builder.num_nodes - nodes_before - num_ancestors_in_epoch)
        mean_memory = np.mean([matcher.total_memory for matcher in self.matcher])
        logger.debug(
            "Finished epoch {} with {} ancestors; {} extra nodes inserted; "
            "mean_tb_size={:.2f} edges={}; mean_matcher_mem={}".format(
                current_time, num_ancestors_in_epoch, extra_nodes,
                np.sum(self.mean_traceback_size) / np.sum(self.num_matches),
                self.tree_sequence_builder.num_edges,
                humanize.naturalsize(mean_memory, binary=True)))
        self.mean_traceback_size[:] = 0
        self.num_matches[:] = 0
        self.results.clear()

    def __match_ancestors_single_threaded(self):
        for j in range(self.start_epoch, self.num_epochs):
            self.__update_progress_epoch(j)
            start, end = map(int, self.epoch_slices[j])
            for ancestor_id in range(start, end):
                a = next(self.ancestors)
                self.__ancestor_find_path(ancestor_id, a)
            self.__complete_epoch(j)

    def __match_ancestors_multi_threaded(self, start_epoch=1):
        # See note on match samples multithreaded below. Should combine these
        # into a single function. Possibly when trying to make the thread
        # error handling more robust.

        queue_depth = 8 * self.num_threads  # Seems like a reasonable limit
        match_queue = queue.Queue(queue_depth)

        def match_worker(thread_index):
            while True:
                work = match_queue.get()
                if work is None:
                    break
                ancestor_id, a = work
                self.__ancestor_find_path(ancestor_id, a, thread_index)
                match_queue.task_done()
            match_queue.task_done()

        match_threads = [
            threads.queue_consumer_thread(
                match_worker, match_queue, name="match-worker-{}".format(j),
                index=j)
            for j in range(self.num_threads)]
        logger.info("Started {} match worker threads".format(self.num_threads))

        for j in range(self.start_epoch, self.num_epochs):
            self.__update_progress_epoch(j)
            start, end = map(int, self.epoch_slices[j])
            for ancestor_id in range(start, end):
                a = next(self.ancestors)
                match_queue.put((ancestor_id, a))
            # Block until all matches have completed.
            match_queue.join()
            self.__complete_epoch(j)

        # Stop the the worker threads.
        for j in range(self.num_threads):
            match_queue.put(None)
        for j in range(self.num_threads):
            match_threads[j].join()

    def match_ancestors(self):
        logger.info("Starting ancestor matching for {} epochs".format(self.num_epochs))
        if self.num_threads <= 0:
            self.__match_ancestors_single_threaded()
        else:
            self.__match_ancestors_multi_threaded()
        ts = self.store_output()
        logger.info("Finished ancestor matching")
        return ts

    def store_output(self):
        if self.num_ancestors > 0:
            ts = self.get_tree_sequence(rescale_positions=False)
        else:
            # Allocate an empty tree sequence.
            ts = msprime.load_tables(
                nodes=msprime.NodeTable(), edges=msprime.EdgeTable(), sequence_length=1)
        if self.output_path is not None:
            ts.dump(self.output_path)
        return ts


class SampleMatcher(Matcher):
    progress_bar_description = "match-samples"

    def __init__(self, sample_data, ancestors_ts, **kwargs):
        super().__init__(sample_data, **kwargs)
        self.restore_tree_sequence_builder(ancestors_ts)
        self.sample_haplotypes = self.sample_data.haplotypes()
        self.sample_ids = np.zeros(self.num_samples, dtype=np.int32)
        for j in range(self.num_samples):
            self.sample_ids[j] = self.tree_sequence_builder.add_node(0)
        self.allocate_progress_monitor(self.num_samples)

    def __process_sample(self, sample_id, haplotype, thread_index=0):
        # print("process sample", haplotype)
        # print("mutated_sites = ", self.mutated_sites)
        # mask = np.zeros(self.num_sites, dtype=np.uint8)
        # mask[self.mutated_sites] = 1
        # h = np.logical_and(haplotype, mask).astype(np.uint8)
        # diffs = np.where(h != haplotype)[0]
        self._find_path(sample_id, haplotype, 0, self.num_sites, thread_index)
        match = self.match[thread_index]
        diffs = np.where(haplotype != match)[0]
        derived_state = haplotype[diffs]
        self.results.set_mutations(sample_id, diffs.astype(np.int32), derived_state)

    def __match_samples_single_threaded(self):
        j = 0
        for a in self.sample_haplotypes:
            sample_id = self.sample_ids[j]
            self.__process_sample(sample_id, a)
            j += 1
        assert j == self.num_samples

    def __match_samples_multi_threaded(self):
        # Note that this function is not almost identical to the match_ancestors
        # multithreaded function above. All we need to do is provide a function
        # to do the matching and some producer for the actual items and we
        # can bring this into a single function.

        queue_depth = 8 * self.num_threads  # Seems like a reasonable limit
        match_queue = queue.Queue(queue_depth)

        def match_worker(thread_index):
            while True:
                work = match_queue.get()
                if work is None:
                    break
                sample_id, a = work
                self.__process_sample(sample_id, a, thread_index)
                match_queue.task_done()
            match_queue.task_done()

        match_threads = [
            threads.queue_consumer_thread(
                match_worker, match_queue, name="match-worker-{}".format(j),
                index=j)
            for j in range(self.num_threads)]
        logger.info("Started {} match worker threads".format(self.num_threads))

        for sample_id, a in zip(self.sample_ids, self.sample_haplotypes):
            match_queue.put((sample_id, a))

        # Stop the the worker threads.
        for j in range(self.num_threads):
            match_queue.put(None)
        for j in range(self.num_threads):
            match_threads[j].join()

    def match_samples(self):
        logger.info("Started matching for {} samples".format(self.num_samples))
        if self.sample_data.num_variant_sites > 0:
            if self.num_threads <= 0:
                self.__match_samples_single_threaded()
            else:
                self.__match_samples_multi_threaded()
            self.progress_monitor.close()
            progress_monitor = tqdm.tqdm(
                desc="update paths", total=self.num_samples, disable=not self.progress)
            for j in range(self.num_samples):
                sample_id = int(self.sample_ids[j])
                left, right, parent = self.results.get_path(sample_id)
                self.tree_sequence_builder.add_path(
                    sample_id, left, right, parent, compress=self.path_compression)
                site, derived_state = self.results.get_mutations(sample_id)
                self.tree_sequence_builder.add_mutations(sample_id, site, derived_state)
                progress_monitor.update()
            progress_monitor.close()
        logger.info("Finished sample matching")

    def finalise(self, simplify=True, stabilise_node_ordering=False):
        logger.info("Finalising tree sequence")
        ts = self.get_tree_sequence(all_sites=True)
        if simplify:
            logger.info("Running simplify on {} nodes and {} edges".format(
                ts.num_nodes, ts.num_edges))
            if stabilise_node_ordering:
                # Ensure all the node times are distinct so that they will have
                # stable IDs after simplifying. This could possibly also be done
                # by reversing the IDs within a time slice. This is used for comparing
                # tree sequences produced by perfect inference.
                tables = ts.tables
                time = tables.nodes.time
                for t in range(1, int(time[0])):
                    index = np.where(time == t)[0]
                    k = index.shape[0]
                    time[index] += np.arange(k)[::-1] / k
                tables.nodes.set_columns(flags=tables.nodes.flags, time=time)
                msprime.sort_tables(**tables.asdict())
                ts = msprime.load_tables(**tables.asdict())
            ts = ts.simplify(
                samples=self.sample_ids, filter_zero_mutation_sites=False)
            logger.info("Finished simplify; now have {} nodes and {} edges".format(
                ts.num_nodes, ts.num_edges))
        return ts


class ResultBuffer(object):
    """
    A wrapper for numpy arrays representing the results of a copying operations.
    """
    def __init__(self):
        self.paths = {}
        self.mutations = {}
        self.lock = threading.Lock()

    def clear(self):
        """
        Clears this result buffer.
        """
        self.paths.clear()
        self.mutations.clear()

    def set_path(self, node_id, left, right, parent):
        with self.lock:
            assert node_id not in self.paths
            self.paths[node_id] = left, right, parent

    def set_mutations(self, node_id, site, derived_state=None):
        if derived_state is None:
            derived_state = np.ones(site.shape[0], dtype=np.uint8)
        with self.lock:
            self.mutations[node_id] = site, derived_state

    def get_path(self, node_id):
        return self.paths[node_id]

    def get_mutations(self, node_id):
        return self.mutations[node_id]
