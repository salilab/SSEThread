"""Make tar.gz files of PDBs for each cluster"""

import tarfile
import os

for cluster in range(2):
    with open('../modeling/cluster%d_ids.dat' % cluster) as fh:
        with tarfile.open('cluster%d_pdbs.tar.gz' % cluster, 'w:gz') as fh_out:
            for n, line in enumerate(fh.readlines()):
                if n % 20 == 0:
                    print("Written %d models" % n)
                pdb = "../modeling/bsms_pdbs/bsms_%d.pdb" % int(line)
                fh_out.add(name=pdb, arcname=os.path.basename(pdb))
