#!/usr/bin/env python

import unittest
import subprocess
import shutil
import os
import sys

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), '..'))

class Tests(unittest.TestCase):
    def assertExists(self, fname):
        self.assertTrue(os.path.exists(fname), msg="%s does not exist" % fname)

    def test_yeast(self):
        """Test model building and analysis of yeast complex"""
        os.chdir(os.path.join(TOPDIR, 'yeast', 'scripts'))
        p = subprocess.check_call(["python", 'sample_multires.py', "--test"])
        for out in ('../outputs/models_0.0.rmf', '../outputs/stat_0.dat'):
            self.assertExists(out)

        p = subprocess.check_call(["python", 'analyze.py'])
        for out in ('score_distribution.pdf', 'violations.pdf',
                    'violation_counts.pdf', 'ContactMap_cluster.pdf',
                    'ccl1.mrc', 'rad3.mrc', 'ssl2.mrc', 'tfb2.mrc',
                    'tfb4.mrc', 'kin28.mrc', 'ssl1.mrc', 'tfb1.mrc',
                    'tfb3.mrc', 'tfb5.mrc'):
            self.assertExists(out)

    def test_human(self):
        """Test model building and analysis of human complex"""
        os.chdir(os.path.join(TOPDIR, 'human', 'scripts'))
        p = subprocess.check_call(["python", 'sample_multires.py', "--test"])
        for out in ('../outputs/models_0.0.rmf', '../outputs/stat_0.dat'):
            self.assertExists(out)

        p = subprocess.check_call(["python", 'analyze.py'])
        for out in ('scores.out', 'scores.pdf',
                    'score_vs_number_viol_restraints.png',
                    'bar_xlink_violation_counts_run1.png',
                    'violation_counts_alltopscoring_34A.txt',
                    'cluster_run1.pkl', 'ContactMap_all_Matrix_CM.pdf',
                    'ccl1.mrc', 'rad3.mrc', 'ssl2.mrc', 'tfb2.mrc',
                    'tfb4.mrc', 'kin28.mrc', 'ssl1.mrc', 'tfb1.mrc',
                    'tfb3.mrc', 'tfb5.mrc'):
            self.assertExists(out)

        p = subprocess.check_call(["python", 'clustering.py'])
        self.assertExists('cluster_run1.out')
        for out in ('cluster_run1.out', 'clusters.pdf'):
            self.assertExists(out)

if __name__ == '__main__':
    unittest.main()
