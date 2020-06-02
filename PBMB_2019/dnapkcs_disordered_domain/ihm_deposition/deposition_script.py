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


# Will finish citation when manuscript is published
'''
system.citations.append(ihm.Citation(
          pmid='25139911', title="SSEThread: Integrative threading of the DNA-PKcs sequence based on data "
          "from chemical cross-linking and hydrogen deuterium exchange",
          journal="Progress in Biophysics and Molecular Biology", 
          year=2019,
          authors=['Saltzberg DJ', 'Hepburn M', 'Pilla KB',
                   'Schriemer D', 'Lees-Miller SP', 'Blundell T', 'Sali A']))
'''



# PSIPRED was used to provide a prediction of the secondary structure for 
# each residue in sequence
system.software.append(ihm.Software(
          name='PSIPRED', classification='secondary structure prediction',
          description='Protein secondary structure prediction based on '
                      'position-specific scoring matrices',
          version='4.0',
          location='http://bioinf.cs.ucl.ac.uk/psipred/'))

# We used various tools from IMP
system.software.append(ihm.Software(
          name="Integrative Modeling Platform (IMP)",
          version="2.2",
          classification="integrative model building",
          description="integrative model building",
          location='https://integrativemodeling.org'))

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
  def add_atoms(self, pdb_file):
    p = Bio.PDB.PDBParser()
    s = p.get_structure('rep', pdb_file)
    for model in s:
        for nchain, chain in enumerate(model):
            #print(nchain, chain)
            for residue in chain.get_residues():
                for atom in residue:
                    #print(residue, atom)
                    coord = atom.get_vector()
                    self._atoms.append(
                                      ihm.model.Atom(asym_unit=asym, 
                                      seq_id = residue.get_id()[1], 
                                      atom_id = "CA", 
                                      type_symbol = "C",
                                      x = coord[0], 
                                      y = coord[1], 
                                      z = coord[2],)
                                      )
    return self._atoms

  def get_atoms(self):
    return self._atoms

start_model = StartingModel(asym, pdb_dataset, 'A')
atoms = start_model.add_atoms("5luq_A_CA.pdb")


rep = ihm.representation.Representation(
        [ihm.representation.ResidueSegment(asym, rigid=True, primitive="sphere", starting_model=start_model)])

# Set up crosslinking restraints. First, we define new cross-linking reagents
bsp = ihm.ChemDescriptor('BSP', chemical_name='Bis(succinimidyl) penta(ethylene glycol)',
              smiles='C1C(N(C(C1)=O)OC(CCOCCOCCOCCOCCOCCC(=O)ON2C(CCC2=O)=O)=O)=O',
              inchi='InChI=1S/C22H32N2O13/c25-17-1-2-18(26)23(17)36-21(29)5-7-31-9-11-33-13-15-35-16-14-34-12-10-32-8-6-22(30)37-24-19(27)3-4-20(24)28/h1-16H2',
              inchi_key='FTYUGLBWKRXQBD-UHFFFAOYSA-N')

# Second, we add the cross-linking files from MSStudio output
l_bsp = ihm.location.InputFileLocation("./BSP-350X-SEC.csv", details="BSP Crosslink File")
l_dss = ihm.location.InputFileLocation("./DSS-350X-SEC.csv", details="DSS Crosslink File")
l_dsg = ihm.location.InputFileLocation("./DSG-Crosslinks-2.csv", details="DSG Crosslink File")

bsp_dataset = ihm.dataset.CXMSDataset(l_bsp)
dss_dataset = ihm.dataset.CXMSDataset(l_dss)
dsg_dataset = ihm.dataset.CXMSDataset(l_dsg)

dsg_restraint = ihm.restraint.CrossLinkRestraint(dsg_dataset, ihm.cross_linkers.dsg)
dss_restraint = ihm.restraint.CrossLinkRestraint(dss_dataset, ihm.cross_linkers.dss)
bsp_restraint = ihm.restraint.CrossLinkRestraint(bsp_dataset, bsp)
system.restraints.extend((dsg_restraint, dss_restraint, bsp_restraint))

xlrs = [dsg_restraint, dss_restraint, bsp_restraint]
xlr_dists = [42, 32, 27]

# Create IHMFeatureList
# - entity_type 'polymer'
# - feture_type 'residue'
#fl = ihm.restraint.FeatureList(entity_type = 'polymer', feature_type='residue')

#rg = ihm.restraint.RestraintGroup()

#system.restraints.extend(rg)


# Add the individual crosslinks to each restraint
f_num=0
for r in range(len(xlrs)):
    for line in open(xlrs[r].dataset.location.path, "r"):
        if "Decision" in line:
            ix = line.split(",").index("Selected Sites")
        if "Accepted" in line:
            ssites = line.split(",")[ix].split(";")[0].replace("[", "").replace("]", "")
            res1 = int(ssites.split("-")[0][1:])
            res2 = int(ssites.split("-")[1][1:])

            distance = ihm.restraint.UpperBoundDistanceRestraint(xlr_dists[r])

            # We only utilized crosslinks with endpoints in the 199-residue disordered region (2575-2774)
            if (res1 > 2574 and res1 <2775) or (res2 > 2574 and res2 <2775):
                # First, create the poly_residue_feature and pseudo_site for each residue
                # add them to the feature
                id1 = f_num
                id2 = f_num+1

                prf1 = ihm.restraint.ResidueFeature(ranges=[asym(res1,res1)])
                prf2 = ihm.restraint.ResidueFeature(ranges=[asym(res2,res2)])
                prf1_id = id1
                prf2_id = id2

                # Make a function to get this pseudosite
                ps1c, ps2c = get_pseudo_site_endpoints(res1, res2, atoms)
                ps1 = ihm.restraint.PseudoSiteFeature(x=ps1c[0], y=ps1c[1], z=ps1c[2], radius=0.0)
                ps2 = ihm.restraint.PseudoSiteFeature(x=ps2c[0], y=ps2c[1], z=ps2c[2], radius=0.0)
                ps1.set_details('This pseudo site corresponds to residue '+str(res1)+' from entity 1. It is not part of the model representation but is part of the experimental restraints')
                ps2.set_details('This pseudo site corresponds to residue '+str(res2)+' from entity 1. It is not part of the model representation but is part of the experimental restraints')
                ps1.set_entity_type("polymer")
                ps2.set_entity_type("polymer")
                ps1._id = id1
                ps2._id = id2
                f_num+=2

                '''
                f1 = ihm.restraint.Feature()
                f2 = ihm.restraint.Feature()
                f1._id = id1
                f2._id = id2
                f1.type = 'pseudo site'
                f1.type = 'pseudo site'
                f1._get_entity_type = 'polymer'
                f2._get_entity_type = 'polymer'
                f1.details = 'This pseudo site corresponds to residue '+str(res1)+' from entity 1. It is not part of the model representation but is part of the experimental restraints'
                f2.details = 'This pseudo site corresponds to residue '+str(res2)+' from entity 1. It is not part of the model representation but is part of the experimental restraints'
                '''

                system.orphan_features.append(prf1)
                system.orphan_features.append(prf2)
                system.orphan_features.append(ps1)
                system.orphan_features.append(ps2)
                # Second, create the derived distance restraints using these features
                #expxl = ihm.restraint.ExperimentalCrossLink(entity.residue(res1), entity.residue(res2))
                #xls.append(expxl)
                ddrxl = ihm.restraint.DerivedDistanceRestraint(dataset=xlrs[r].dataset, 
                                                            feature1 = ps1, 
                                                            feature2 = ps2, 
                                                            distance = distance)
                system.restraints.append(ddrxl)
                #rxl = ihm.restraint.ResidueCrossLink(expxl, asym, asym, distance)
                #xlxls.append(rxl)
                #xlrs[r].cross_links.append(rxl)
                #xlrs[r].experimental_cross_links.append([expxl])
                #xlrs[r].experimental_cross_links.append([ddrxl])
                #xlrs[r].cross_links.append(rxl)

# Now we add information about how the modeling was done by defining one
# or more protocols. Here we enumerated all possible threading solutions
# and chose the top 10000 scoring models
all_datasets = ihm.dataset.DatasetGroup((dsg_dataset,
                                         dss_dataset, bsp_dataset))

protocol = ihm.protocol.Protocol(name='Modeling')
protocol.steps.append(ihm.protocol.Step(
                        assembly=assembly,
                        dataset_group=all_datasets,
                        method='Enumeration',
                        name='Production sampling',
                        num_models_begin=0,
                        num_models_end=1000, multi_scale=False))

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



# Here, we place a subset of the best scoring models into the deposition
models = []
for pdb_file in glob.glob("../modeling/bsms_pdbs/bsms_*100.pdb")[0:1]:
    #print(pdb_file)
    m = Model(assembly = assembly, protocol=protocol, representation=rep,
        file_name=pdb_file, asym_units=[asym])
    m.get_residues()
    models.append(m)

#---------------
# Write .dcd file of all of the best scoring models
'''
dcd = 'best_scoring_models.dcd'
fh = open(dcd, "wb")
d = ihm.model.DCDWriter(fh)

for pdb_file in glob.glob("../modeling/bsms_pdbs/bsms_*100.pdb"):
    m = Model(assembly = assembly, protocol=protocol, representation=rep,
        file_name=pdb_file, asym_units=[asym])
    d.add_model(m)
'''
l_dcd = ihm.location.OutputFileLocation('best_scoring_models.dcd')

mg = ihm.model.ModelGroup(models)
me = ihm.model.Ensemble(mg, len(models),
    post_process="filtering",
    file=l_dcd)

# Groups are then placed into states, which can in turn be grouped. In this
# case we have only a single state:
state = ihm.model.State([mg],
                    type='Threading Ensemble',
                    name='Threading Ensemble Solution')
system.state_groups.append(ihm.model.StateGroup([state]))

# Once the system is complete, we can write it out to an mmCIF file:
with open('output.cif', 'w') as fh:
    ihm.dumper.write(fh, [system])
