## Heap based Segment Tree optimized for Stab Queries
#
# Author: Tanner Robart
# Date:	  October, 2019

"""
Problem:

The output of an NGS sequencer is a large collection of such disjointed sequences, called
"reads". For the purpose of this puzzle, each read has two properties: (1) a start position within the genome and (2) its
length. Both properties are expressed in base pair units, i.e. letters of DNA like A, T ,C,G. For example, the read "AAATCGA"
has length 7.
A key component of making sense of NGS data is calculating "coverage" (also known as "read depth") at a given genomic
position. At its most basic, coverage is simply the number of reads overlapping a position in the genome. High coverage
gives a measure of confidence in the sequencing results, and the calculation of coverage is basic to sequencing.

The exercise is to write a computer program in Python3 (only using the standard library) to calculate the coverage
for each of the positions in loci.csv, based on the reads in reads.csv.

Approach:

This problem seems to be a stabbing query specific subproblem of the interval and segment tree problems.
Specifically, we want to minimize the number of reads we check in order to count the
number of overlapping segments given a integer position. Using a segment tree
(binary tree that stores segment information in its inherent structure)
should allow us to count the number of overlapping intervals in O(k + log(n)) time,
where n is the number of distinct intervals (smaller than our N segments) and k is the number of retrieved intervals during a query,
in this case, the height of our tree which is log(n), giving us O(2log(n)) ~>  O(log(n))

However, by using a heap structure and only storing the counts in each node,
instead of the segments themselves, we can store the values in an arbitrary data structue
for later access based on its node index. we will use an array.
By doing this, we can reduce the storage cost of the tree significantly, and stab query complexity cost down to O(log(n)),
The optimal possible solution.
"""

import random
import bisect
import os
import math
import csv
import unittest
import pickle
import timeit
from tqdm import tqdm  # non-base python but just for nice progress bars!

input_file = "reads.csv"
output_file = "loci.csv"
points_file = "distinct_points.txt"
# pickleTree = "bob.pkl"
pickleTree = "segmentTree.pkl"


class SegmentTree:
    ## Heap based segment tree implemented as an array
    def __init__(self, points, segments, node=1, uncompressed=False, _n=0):
        n = len(points)
        if _n > 0:
            n = _n
        self.__total = 2 * n - 1
        self.__tree = [0] * int(self.__total + 2)  # actual tree data array of counts
        self.__N = int(self.__total / 2) - 1
        self.__points = points

        # Insertion of the intervals into the tree takes O(n * log(n)) time

        for i in tqdm(segments, desc="inserting segments:"):
            startpoint, endpoint = i
            if not uncompressed:
                startpoint, endpoint = (
                    bisect.bisect_right(self.__points, startpoint) + 1,
                    bisect.bisect_right(self.__points, endpoint) + 1,
                )
            self.insert(startpoint, endpoint)

    def insert(self, startpoint, endpoint):
        # Bottom up way to increment count for every node from the leaves
        # up to the elementary interval node (the node at which the
        # union of the two start and end leaf nodes occurs)
        # complexity O(log(n))
        l = startpoint + self.__N
        r = endpoint + self.__N
        # print("startpoint:", startpoint, "leaf node:", l)
        # print("startpoint:", startpoint, "leaf node:", r)
        while l <= r:
            if l % 2 == 1:
                self.__tree[l] += 1
            l = int(math.floor((l + 1) / 2))
            if r % 2 == 0:
                self.__tree[r] += 1
            r = int(math.floor((r - 1) / 2))

    def stab_query(self, point):
        # Bottom up way to return the sum of counts stored in
        # one traversal from leaf to root of the tree.
        # complexity O(log(n))
        count = 0

        # perform a binary search of the ordered array to find
        # the index of the leaf node corresponding to the given point
        node = int(bisect.bisect_right(self.__points, point) + self.__N) + 1

        while not self.complete(node):
            count += self.__tree[node]
            node = int(math.floor(node / 2))
            if node == 0:
                break
        return count

    def stab_query_uncompressed(self, point):  # gausa binser, gausah ada array bahkan
        # Bottom up way to return the sum of counts stored in
        # one traversal from leaf to root of the tree.
        # complexity O(log(n))
        count = 0
        node = point + 1 * self.__N
        while not self.complete(node):
            count += self.__tree[node]
            node = int(math.floor(node / 2))
            if node == 0:
                break
        return count

    def complete(self, node):
        # Boolean check if node is a complete binary subtree
        formula = 2 ** (
            math.floor(math.log(self.lowest_incomplete(node)))
            - math.floor(math.log(node))
        )
        z = math.floor(self.lowest_incomplete(node) / formula)
        if node == z:
            return True
        else:
            return False

    def lowest_incomplete(self, node):
        # Returns the lowest incomplete node
        t = math.log(self.__N & -self.__N)
        return int(math.floor(float(self.__N) / (2 ** (1 + t))))


def load_points(points_file):

    # Preprocessing: loads all distinct points from reads.csv,
    # then sorts all points, the sort takes O(n*log(n)) time using merge sort,
    # where n is the number of distinct start and end points.
    # Does this once and saves the list of points for use, to avoid re-sorting a 4~ish million item list again
    if not os.path.isfile(points_file):
        with open(input_file, "r") as csvfile:
            points = set()  # using set to reduce to distinct values
            reader = csv.DictReader(csvfile)
            for row in tqdm(reader, desc="building points list"):
                startpoint, length = int(row["start"]), int(row["length"])
                points.add(startpoint)
                points.add(startpoint + length)
        points = list(points)
        points.sort()
        with open(points_file, "w") as f:
            f.writelines("%s\n" % l for l in points)
    else:
        # otherwise load presorted distinct points list
        points = [int(line.rstrip("\n")) for line in open(points_file)]

    return points


def main():

    # check if we have already constructed and pickled the segment tree
    if not os.path.isfile(pickleTree):

        # performs preproccessing of points, or loads if already done
        points = load_points(points_file)

        # open reads.csv file and extract segments
        with open(input_file, "r") as csvfile:
            segments = []
            reader = csv.DictReader(csvfile)
            for row in reader:
                # use length to convert into start and end points all in one big array
                # for preprocessing for effienciently sorting and building segment tree,
                # and for storing segment values
                startpoint, length = int(row["start"]), int(row["length"])
                segments.append((startpoint, startpoint + length))

        # intiialize a balanced binary search tree with the number of leaves are
        # equal to distinct start and endpoints (in 'points' object), then populate it by inserting all the segments
        print(points, segments)
        segmentTree = SegmentTree(points, segments)
        print(segmentTree._SegmentTree__tree)

        # free up the segments and points objects are
        # now implicitly encoded within our segment tree
        segments = None
        points = None

        with open(pickleTree, "wb") as f:
            pickle.dump(segmentTree, f)
            print("Data has been pickled to file {}".format(f))

    else:
        # pickle load the already constructed segment tree
        with open(pickleTree, "rb") as f:
            print("Loading pickled segment tree, yum!")
            segmentTree = pickle.load(f)
            # print(segmentTree._SegmentTree__tree)

    # Open loci.csv to get positions for stab queries,
    # make stab queries and store the returned interval counts
    with open(output_file, "r") as csvfile:
        loci_list = [("position", "coverage")]
        reader = csv.DictReader(csvfile)
        for row in reader:
            position = int(row["position"])
            coverage = segmentTree.stab_query(position)
            loci_list.append((position, coverage))

    # Write interval counts to loci.csv
    with open(output_file, "w") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(loci_list)
        print("wrote to loci.csv!")
    return loci_list


class SolutionTest(unittest.TestCase):
    def setUp(self):
        pass

    # @unittest.skip("kalau mau test yg uncompressed aja")
    def test_toy_set(self):
        # toy test set for debugging
        # points = random.sample(range(10000, 1000000), 100000)
        # segments = random.sample(points, 9000)

        # segments = [(i, random.sample(points, 1)[0]) for i in segments]
        # segments = [i for i in segments if i[0] > i[1]]
        segments = []
        segments.extend(
            [
                (3, 4),
                (3, 4),
                (3, 4),
                (3, 4),
                (3, 6),
                (10, 20),
                (6, 50),
                (3, 3),
                (1000, 7000000),
            ]
        )
        points = [3, 4, 6, 7, 8, 9, 10, 20, 50, 76, 89, 1000, 7000000]
        points.sort()

        segmentTree = SegmentTree(points, segments)
        self.assertEqual(2, segmentTree.stab_query(15), "output for 15 should be 2")
        self.assertEqual(6, segmentTree.stab_query(3), "output for 3 should be: 6")
        self.assertEqual(5, segmentTree.stab_query(4), "output for 4 should be: 5")
        self.assertEqual(1, segmentTree.stab_query(7), "output for 7 should be: 1")
        self.assertEqual(2, segmentTree.stab_query(6), "output for 6 should be: 2")
        self.assertEqual(
            1, segmentTree.stab_query(10000), "output for 10000 should be: 1"
        )

    @unittest.skip("kalau mau test yg compressed aja")
    def test_toy_set_uncompressed(self):
        segments = []
        segments.extend(
            [
                (3, 4),
                (3, 4),
                (3, 4),
                (3, 4),
                (3, 6),
                (10, 20),
                (6, 50),
                (3, 3),
                (1000, 7000000),
            ]
        )
        points = [3, 4, 6, 7, 8, 9, 10, 20, 50, 76, 89, 1000, 7000000]
        # points = [i for i in range(1, max(points) + 1)] # ga perlu array lagi
        # points.sort() # gaperlu sort lagi

        segmentTree = SegmentTree(points, segments, 1, 1, max(points))
        self.assertEqual(
            2, segmentTree.stab_query_uncompressed(15), "output for 15 should be 2"
        )
        self.assertEqual(
            6, segmentTree.stab_query_uncompressed(3), "output for 3 should be: 6"
        )
        self.assertEqual(
            5, segmentTree.stab_query_uncompressed(4), "output for 4 should be: 5"
        )
        self.assertEqual(
            1, segmentTree.stab_query_uncompressed(7), "output for 7 should be: 1"
        )
        self.assertEqual(
            2, segmentTree.stab_query_uncompressed(6), "output for 6 should be: 2"
        )
        self.assertEqual(
            1,
            segmentTree.stab_query_uncompressed(10000),
            "output for 10000 should be: 1",
        )


if __name__ == "__main__":

    unittest.main()
