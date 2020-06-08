#!/usr/bin/env python

import unittest
import os
import shutil
import sys
import subprocess
import ihm.reader

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), '..'))

class Tests(unittest.TestCase):
    def test_mmcif(self):
        """Test generation of mmCIF output"""
        os.chdir(os.path.join(TOPDIR, 'PBMB_2019', 'dnapkcs_disordered_domain',
                              'ihm_deposition'))
        if os.path.exists("output.cif"):
            os.unlink("output.cif")
        # Potentially override methods that need network access
        env = os.environ.copy()
        env['PYTHONPATH'] = os.path.join(TOPDIR, 'test', 'mock') \
                            + ':' + env.get('PYTHONPATH', '')
        p = subprocess.check_call(
                [sys.executable, "deposition_script.py"], env=env)
        # Check output file
        self._check_mmcif_file('output.cif')

    def _check_mmcif_file(self, fname):
        with open(fname) as fh:
            s, = ihm.reader.read(fh)
        self.assertEqual(len(s.citations), 1)
        self.assertEqual(s.citations[0].doi, '10.1016/j.pbiomolbio.2019.09.003')
        self.assertEqual(len(s.software), 3)
        self.assertEqual(len(s.orphan_starting_models), 1)
        # Should be 1 state
        self.assertEqual(len(s.state_groups), 1)
        state1, = s.state_groups[0]
        # Should be 10 models
        self.assertEqual(sum(len(x) for x in state1), 10)
        # Check # of spheres and atoms in each model
        for m in state1[0]:
            self.assertEqual(len(m._spheres), 3600)
            self.assertEqual(len(m._atoms), 0)
        # Should be 2 ensembles
        self.assertEqual([e.num_models for e in s.ensembles], [2758, 2242])
        # 3 different crosslink restraints
        xl1, xl2, xl3 = s.restraints
        self.assertEqual(len(xl1.experimental_cross_links), 63)
        self.assertEqual(len(xl1.cross_links), 20)
        self.assertEqual(xl1.linker.auth_name, 'DSG')

        self.assertEqual(len(xl2.experimental_cross_links), 48)
        self.assertEqual(len(xl2.cross_links), 12)
        self.assertEqual(xl2.linker.auth_name, 'DSS')

        self.assertEqual(len(xl3.experimental_cross_links), 52)
        self.assertEqual(len(xl3.cross_links), 14)
        self.assertEqual(xl3.linker.auth_name, 'BSP')


if __name__ == '__main__':
    unittest.main()
