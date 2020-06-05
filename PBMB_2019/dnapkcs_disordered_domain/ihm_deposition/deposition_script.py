# This example demonstrates the use of the Python IHM library to generate
# an mmCIF file for a very simple integrative docking study. Two subunits,
# A and B, each of which is fitted against small angle X-ray (SAXS) data, are
# docked together into a complex, AB, which is fitted against an electron
# microscopy density map.
import sys
sys.path.append("./python-ihm")
import ihm
import ihm.location
import ihm.dataset
import ihm.representation
import ihm.restraint
import ihm.protocol
import ihm.analysis
import ihm.model
import ihm.dumper
import ihm.startmodel
import ihm.cross_linkers
import Bio.PDB
import Bio.SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
import glob
from DisorderedCrossLink import *

# First, we create a system, which contains everything we know about the
# modeling. A single mmCIF file can contain multiple Systems, but in most
# cases we use just one:

title = ("Integrative threading of the DNA-PKcs sequence based on data from "
  "chemical cross-linking and hydrogen deuterium exchange")
system = ihm.System(title=title)

system.citations.append(ihm.Citation.from_pubmed_id(31570166))

# PSIPRED was used to provide a prediction of the secondary structure for 
# each residue in sequence
system.software.append(ihm.Software(
          name='PSIPRED', classification='secondary structure prediction',
          description='Protein secondary structure prediction based on '
                      'position-specific scoring matrices',
          version='4.0',
          location='http://bioinf.cs.ucl.ac.uk/psipred/'))

# We used various tools from IMP
imp_software = ihm.Software(
          name="Integrative Modeling Platform (IMP)",
          version="2.2",
          classification="integrative model building",
          description="integrative model building",
          location='https://integrativemodeling.org')
system.software.append(imp_software)

# We used scikit-learn for clustering
sklearn_software = ihm.Software(
          name="scikit-learn",
          version="0.21.3",
          classification="model building",
          description="Machine learning in Python",
          location='https://scikit-learn.org/stable/')
system.software.append(sklearn_software)

# Next, we describe the input data we used, using dataset classes.
pdb_l = ihm.location.PDBLocation("5LUQ", version="1.1")
pdb_dataset = ihm.dataset.PDBDataset(pdb_l)

# Next, define the entities for each unique sequence in the system
for record in Bio.SeqIO.parse("P78527.fasta", "fasta"):
    sequence = record.seq

entity = ihm.Entity(sequence, description='DNA-PKcs')
system.entities.append(entity)

# Next, we define asymmetric units for everything we modeled.
asym = ihm.AsymUnit(entity, details='Threading Model')
system.asym_units.append(asym)

# Next, we group asymmetric units (and/or entities) into assemblies.
assembly = ihm.Assembly([asym], name='Threading Model')


# The system was represented as a bead model with one residue per bead
# Using chain A of 5LUQ pdb
class StartingModel(ihm.startmodel.StartingModel):
    def __init__(self, pdb_file, **kwargs):
        super(StartingModel, self).__init__(**kwargs)
        self.pdb_file = pdb_file

    def get_atoms(self):
        p = Bio.PDB.PDBParser()
        s = p.get_structure('rep', self.pdb_file)
        for model in s:
            for nchain, chain in enumerate(model):
                #print(nchain, chain)
                for residue in chain.get_residues():
                    for atom in residue:
                        #print(residue, atom)
                        coord = atom.get_vector()
                        yield ihm.model.Atom(
                                asym_unit=asym, seq_id=residue.get_id()[1],
                                atom_id="CA", type_symbol="C",
                                x=coord[0], y=coord[1], z=coord[2])


start_model = StartingModel(asym_unit=asym, dataset=pdb_dataset, asym_id='A',
                            pdb_file="../data/5luq_A_CA.pdb")


rep = ihm.representation.Representation(
        [ihm.representation.ResidueSegment(asym, rigid=True, primitive="sphere", starting_model=start_model)])

# Set up crosslinking restraints. First, we define new cross-linking reagents
bsp = ihm.ChemDescriptor('BSP', chemical_name='Bis(succinimidyl) penta(ethylene glycol)',
              smiles='C1C(N(C(C1)=O)OC(CCOCCOCCOCCOCCOCCC(=O)ON2C(CCC2=O)=O)=O)=O',
              inchi='InChI=1S/C22H32N2O13/c25-17-1-2-18(26)23(17)36-21(29)5-7-31-9-11-33-13-15-35-16-14-34-12-10-32-8-6-22(30)37-24-19(27)3-4-20(24)28/h1-16H2',
              inchi_key='FTYUGLBWKRXQBD-UHFFFAOYSA-N')

# Second, we add the cross-linking files from MSStudio output
l_bsp = ihm.location.InputFileLocation(
        "../data/BSP-350X-SEC.csv", details="BSP Crosslink File")
l_dss = ihm.location.InputFileLocation(
        "../data/DSS-350X-SEC.csv", details="DSS Crosslink File")
l_dsg = ihm.location.InputFileLocation(
        "../data/DSG-Crosslinks-2.csv", details="DSG Crosslink File")

bsp_dataset = ihm.dataset.CXMSDataset(l_bsp)
dss_dataset = ihm.dataset.CXMSDataset(l_dss)
dsg_dataset = ihm.dataset.CXMSDataset(l_dsg)

dsg_restraint = ihm.restraint.CrossLinkRestraint(dsg_dataset, ihm.cross_linkers.dsg)
dss_restraint = ihm.restraint.CrossLinkRestraint(dss_dataset, ihm.cross_linkers.dss)
bsp_restraint = ihm.restraint.CrossLinkRestraint(bsp_dataset, bsp)
system.restraints.extend((dsg_restraint, dss_restraint, bsp_restraint))

xlrs = [dsg_restraint, dss_restraint, bsp_restraint]
xlr_dists = [42, 32, 27]

# Add all experimentally-determined cross-links
for r, dist in zip(xlrs, xlr_dists):
    with open(r.dataset.location.path) as fh:
        distance = ihm.restraint.UpperBoundDistanceRestraint(dist)
        for line in fh:
            if "Decision" in line:
                ix = line.split(",").index("Selected Sites")
            if "Accepted" in line:
                ssites = line.split(",")[ix].split(";")[0].replace("[", "").replace("]", "")
                res1 = int(ssites.split("-")[0][1:])
                res2 = int(ssites.split("-")[1][1:])
                ex_xl = ihm.restraint.ExperimentalCrossLink(
                            entity.residue(res1), entity.residue(res2))
                r.experimental_cross_links.append([ex_xl])
                # We only utilized crosslinks with endpoints in the
                # 199-residue disordered region (2575-2774)
                if ((res1 > 2574 and res1 < 2775)
                    or (res2 > 2574 and res2 < 2775)):
                    # Add cross-link between residues
                    rcl = ihm.restraint.ResidueCrossLink(
                            ex_xl, asym1=asym, asym2=asym, distance=distance)
                    r.cross_links.append(rcl)


def add_model_cross_links(m):
    """Add cross-link information for a given model (the endpoint
       coordinates are model-dependent)."""

    psf = PseudoSiteFinder(m)
    for r in xlrs:
        for xl in r.cross_links:
            ex_xl = xl.experimental_cross_link
            res1 = ex_xl.residue1.seq_id
            res2 = ex_xl.residue2.seq_id
            # Make a function to get this pseudosite
            ps1, ps2, form = psf.get_pseudo_site_endpoints(res1, res2)

            # Add pseudo-site to the existing restraint
            if xl.pseudo1 is None:
                xl.pseudo1 = []
            if xl.pseudo2 is None:
                xl.pseudo2 = []
            xl.pseudo1.append(ps1)
            xl.pseudo2.append(ps2)

# Now we add information about how the modeling was done by defining one
# or more protocols. Here we enumerated all possible threading solutions
# and chose the top 5000 scoring models
all_datasets = ihm.dataset.DatasetGroup((dsg_dataset,
                                         dss_dataset, bsp_dataset))

protocol = ihm.protocol.Protocol(name='Modeling')
modeling_script = ihm.location.WorkflowFileLocation(
        "../modeling/modeling_enumerate_multi.py",
        details="Main modeling script")
protocol.steps.append(ihm.protocol.Step(
                        assembly=assembly,
                        dataset_group=all_datasets,
                        method='Enumeration',
                        name='Production sampling',
                        num_models_begin=0,
                        num_models_end=2860000, multi_scale=False,
                        script_file=modeling_script,
                        software=imp_software))
analysis = ihm.analysis.Analysis()
analysis.steps.append(ihm.analysis.FilterStep(
    feature='energy/score', num_models_begin=2860000, num_models_end=5000,
    assembly=assembly,
    details="Selection of the top 5000 best-scoring solutions"))
cluster_script = ihm.location.WorkflowFileLocation(
        "../modeling/cluster_bsms.py",
        details="Clustering script using KMeans from scikit-learn")
analysis.steps.append(ihm.analysis.ClusterStep(
    feature='RMSD', num_models_begin=5000, num_models_end=5000,
    assembly=assembly, details="Clustering using KMeans",
    script_file=cluster_script, software=sklearn_software))
protocol.analyses.append(analysis)

# Rather than storing another copy of the coordinates in the IHM library
# (which could use a lot of memory), we need to provide a mechanism to
# translate them into the IHM data model. We do this straighforwardly by
# subclassing the IHM Model class and overriding the get_atoms
# and get_spheres methods:

class Model(ihm.model.Model):
    """Pass a BioPython model through to IHM"""
    def __init__(self, file_name, asym_units, **kwargs):
        super(Model, self).__init__(**kwargs)
        self.file_name = file_name
        self.asym_units = asym_units
        self.structured_residues = []

    def get_residues(self):
        # Use BioPython to read the structure from a PDB file, and then yield
        # a set of ihm.model.Residue objects
        p = Bio.PDB.PDBParser()
        s = p.get_structure('rep', self.file_name)

        for model in s:
            for nchain, chain in enumerate(model):
                asym = self.asym_units[nchain]
                for residue in chain.get_residues():
                    for atom in residue:
                        self.structured_residues.append(residue.get_id()[1])
                        #print(residue, atom)
                        coord = atom.get_vector()
                        self.add_sphere(ihm.model.Sphere(asym_unit=asym, seq_id_range=(residue.get_id()[1], residue.get_id()[1]),
                                x=coord[0], y=coord[1],
                                z=coord[2], radius=1.0, rmsf=0.0))


    def get_unmodeled_residue_coordinates(self, ps_res):
      # ps_res = Residue to create the pseudo-site for
      #
      #
      # ep1 = Residue of nearest N-terminal endpoint
      # ep2 = Residue of nearest C-terminal endpoint
      # From the threading model, calculate the position of C'v for all unmodeled residues
      # Returns the xyz for each endpoint plus the 
      xyz1 = ep1[0].get_vector()
      xyz2 = ep2[0].get_vector()

      cv_prime = None
      return coords

    def create_pseudo_site(self, residue):
      # Given a residue, find its position according to the model and 
      # return an ihm.restraint.CrossLinkPseudoSite
      coords = [0,0,0]
      ps = ihm.restraint.PseudoSite(coords[0], coords[1], coords[2])
      xlps = ihm.restraint.CrossLinkPseudoSite(ps, self)
      return xlps

    def create_pseudocrosslink(self, expxl, asym1, asym2, distance, resnum1, resnum2):
      # Each crosslink will have two coordinates. These coordinates are either from a modeled residue
      # in which case, the crosslink coordinates are simple the coordiantes of those residues.
      #
      # For those residues that are unmodeled, we must calculate a pseudosite for the residue. This pseudosite
      # will be different for each crosslink and depends on the position of the nearest modeled residues
      xyz1, xyz2 = get_pseudo_site_endpoints(resnum1, resnum2, atoms)
      xlps1 = ihm.restraint.CrossLinkPseudoSite(ihm.restraint.PseudoSite(xyz1[0], xyz1[1], xyz1[2]), self)
      xlps2 = ihm.restraint.CrossLinkPseudoSite(ihm.restraint.PseudoSite(xyz2[0], xyz2[1], xyz2[2]), self)
      rxl = ihm.restraint.ResidueCrossLink(expxl, asym1, asym2, distance, pseudo1=xlps1, pseudo2=xlps2)
      return rxl 



# Here, we place a subset of the best scoring models for each cluster
# into the deposition
mgs = []
for cluster in range(2):
    with open('../modeling/cluster%d_ids.dat' % cluster) as fh:
        pdbs = ["../modeling/bsms_pdbs/bsms_%d.pdb" % int(x)
                for x in fh.readlines()]
    models = []
    for pdb_file in pdbs[0:1]:
        #print(pdb_file)
        m = Model(assembly = assembly, protocol=protocol, representation=rep,
                  file_name=pdb_file, asym_units=[asym],
                  name="Example model for cluster %d" % cluster)
        m.get_residues()
        add_model_cross_links(m)
        models.append(m)

    #---------------
    # Write .dcd file of all of the best scoring models if script called
    # with --dcd option
    dcd = 'cluster%d.dcd' % cluster
    if '--dcd' in sys.argv:
        with open(dcd, "wb") as fh:
            d = ihm.model.DCDWriter(fh)

            for i, pdb_file in enumerate(pdbs):
                if i % 20 == 0:
                    print("Added %d of %d models to DCD" % (i, len(pdbs)))
                m = Model(assembly=assembly, protocol=protocol,
                          representation=rep, file_name=pdb_file,
                          asym_units=[asym])
                m.get_residues()
                d.add_model(m)
    l_dcd = ihm.location.OutputFileLocation('cluster%d.dcd' % cluster)

    mg = ihm.model.ModelGroup(models, name="Cluster %d" % cluster)
    mgs.append(mg)
    me = ihm.model.Ensemble(mg, len(pdbs),
        post_process=protocol.analyses[-1],
        file=l_dcd, name="Cluster %d" % cluster)
    system.ensembles.append(me)

# Groups are then placed into states, which can in turn be grouped. In this
# case we have only a single state:
state = ihm.model.State(mgs,
                    type='Threading Ensemble',
                    name='Threading Ensemble Solution')
system.state_groups.append(ihm.model.StateGroup([state]))

# Rewrite local paths to point to Zenodo
r = ihm.location.Repository(
    doi="10.5281/zenodo.2580423",
    url="https://zenodo.org/record/2580424/files/salilab/SSEThread-v0.1.zip",
    top_directory="salilab-SSEThread-bf6b912",
    root="../../..")
system.update_locations_in_repositories([r])

# Once the system is complete, we can write it out to an mmCIF file:
with open('output.cif', 'w') as fh:
    ihm.dumper.write(fh, [system])
