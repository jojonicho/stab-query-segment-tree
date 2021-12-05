elf.assertEqual(
            2, segmentTree.stab_query_uncompressed(6), "output for 6 should be: 2"
        )
        self.assertEqual(
            1,
            segmentTree.stab_query_uncompressed(10000),
            "output for 10000 should be: 1",
        )