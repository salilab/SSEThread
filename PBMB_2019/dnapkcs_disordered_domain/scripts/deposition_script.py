# This example demonstrates the use of the Python IHM library to generate
# an mmCIF file for a very simple integrative docking study. Two subunits,
# A and B, each of which is fitted against small angle X-ray (SAXS) data, are
# docked together into a complex, AB, which is fitted against an electron
# microscopy density map.

import ihm
import ihm.location
import ihm.dataset
import ihm.representation
import ihm.restraint
import ihm.protocol
import ihm.model
import ihm.dumper
import ihm.cross_linkers
import Bio.PDB
import glob

# First, we create a system, which contains everything we know about the
# modeling. A single mmCIF file can contain multiple Systems, but in most
# cases we use just one:

title = ("Integrative threading of the DNA-PKcs sequence based on data from chemical cross-linking and hydrogen deuterium exchange")
system = ihm.System(title=title)

sequence = "MAGSGAGVRCSLLRLQETLSAADRCGAALAGHQLIRGLGQECVLSSSPAVLALQTSLVFSRDFGLLVFVRKSLNSIEFRECREEILKFLCIFLEKMGQKIAPYSVEIKNTCTSVYTKDRAAKCKIPALDLLIKLLQTFRSSRLMDEFKIGELFSKFYGELALKKKIPDTVLEKVYELLGLLGEVHPSEMINNAENLFRAFLGELKTQMTSAVREPKLPVLAGCLKGLSSLLCNFTKSMEEDPQTSREIFNFVLKAIRPQIDLKRYAVPSAGLRLFALHASQFSTCLLDNYVSLFEVLLKWCAHTNVELKKAALSALESFLKQVSNMVAKNAEMHKNKLQYFMEQFYGIIRNVDSNNKELSIAIRGYGLFAGPCKVINAKDVDFMYVELIQRCKQMFLTQTDTGDDRVYQMPSFLQSVASVLLYLDTVPEVYTPVLEHLVVMQIDSFPQYSPKMQLVCCRAIVKVFLALAAKGPVLRNCISTVVHQGLIRICSKPVVLPKGPESESEDHRASGEVRTGKWKVPTYKDYVDLFRHLLSSDQMMDSILADEAFFSVNSSSESLNHLLYDEFVKSVLKIVEKLDLTLEIQTVGEQENGDEAPGVWMIPTSDPAANLHPAKPKDFSAFINLVEFCREILPEKQAEFFEPWVYSFSYELILQSTRLPLISGFYKLLSITVRNAKKIKYFEGVSPKSLKHSPEDPEKYSCFALFVKFGKEVAVKMKQYKDELLASCLTFLLSLPHNIIELDVRAYVPALQMAFKLGLSYTPLAEVGLNALEEWSIYIDRHVMQPYYKDILPCLDGYLKTSALSDETKNNWEVSALSRAAQKGFNKVVLKHLKKTKNLSSNEAISLEEIRIRVVQMLGSLGGQINKNLLTVTSSDEMMKSYVAWDREKRLSFAVPFREMKPVIFLDVFLPRVTELALTASDRQTKVAACELLHSMVMFMLGKATQMPEGGQGAPPMYQLYKRTFPVLLRLACDVDQVTRQLYEPLVMQLIHWFTNNKKFESQDTVALLEAILDGIVDPVDSTLRDFCGRCIREFLKWSIKQITPQQQEKSPVNTKSLFKRLYSLALHPNAFKRLGASLAFNNIYREFREEESLVEQFVFEALVIYMESLALAHADEKSLGTIQQCCDAIDHLCRIIEKKHVSLNKAKKRRLPRGFPPSASLCLLDLVKWLLAHCGRPQTECRHKSIELFYKFVPLLPGNRSPNLWLKDVLKEEGVSFLINTFEGGGCGQPSGILAQPTLLYLRGPFSLQATLCWLDLLLAALECYNTFIGERTVGALQVLGTEAQSSLLKAVAFFLESIAMHDIIAAEKCFGTGAAGNRTSPQEGERYNYSKCTVVVRIMEFTTTLLNTSPEGWKLLKKDLCNTHLMRVLVQTLCEPASIGFNIGDVQVMAHLPDVCVNLMKALKMSPYKDILETHLREKITAQSIEELCAVNLYGPDAQVDRSRLAAVVSACKQLHRAGLLHNILPSQSTDLHHSVGTELLSLVYKGIAPGDERQCLPSLDLSCKQLASGLLELAFAFGGLCERLVSLLLNPAVLSTASLGSSQGSVIHFSHGEYFYSLFSETINTELLKNLDLAVLELMQSSVDNTKMVSAVLNGMLDQSFRERANQKHQGLKLATTILQHWKKCDSWWAKDSPLETKMAVLALLAKILQIDSSVSFNTSHGSFPEVFTTYISLLADTKLDLHLKGQAVTLLPFFTSLTGGSLEELRRVLEQLIVAHFPMQSREFPPGTPRFNNYVDCMKKFLDALELSQSPMLLELMTEVLCREQQHVMEELFQSSFRRIARRGSCVTQVGLLESVYEMFRKDDPRLSFTRQSFVDRSLLTLLWHCSLDALREFFSTIVVDAIDVLKSRFTKLNESTFDTQITKKMGYYKILDVMYSRLPKDDVHAKESKINQVFHGSCITEGNELTKTLIKLCYDAFTENMAGENQLLERRRLYHCAAYNCAISVICCVFNELKFYQGFLFSEKPEKNLLIFENLIDLKRRYNFPVEVEVPMERKKKYIEIRKEAREAANGDSDGPSYMSSLSYLADSTLSEEMSQFDFSTGVQSYSYSSQDPRPATGRFRRREQRDPTVHDDVLELEMDELNRHECMAPLTALVKHMHRSLGPPQGEEDSVPRDLPSWMKFLHGKLGNPIVPLNIRLFLAKLVINTEEVFRPYAKHWLSPLLQLAASENNGGEGIHYMVVEIVATILSWTGLATPTGVPKDEVLANRLLNFLMKHVFHPKRAVFRHNLEIIKTLVECWKDCLSIPYRLIFEKFSGKDPNSKDNSVGIQLLGIVMANDLPPYDPQCGIQSSEYFQALVNNMSFVRYKEVYAAAAEVLGLILRYVMERKNILEESLCELVAKQLKQHQNTMEDKFIVCLNKVTKSFPPLADRFMNAVFFLLPKFHGVLKTLCLEVVLCRVEGMTELYFQLKSKDFVQVMRHRDDERQKVCLDIIYKMMPKLKPVELRELLNPVVEFVSHPSTTCREQMYNILMWIHDNYRDPESETDNDSQEIFKLAKDVLIQGLIDENPGLQLIIRNFWSHETRLPSNTLDRLLALNSLYSPKIEVHFLSLATNFLLEMTSMSPDYPNPMFEHPLSECEFQEYTIDSDWRFRSTVLTPMFVETQASQGTLQTRTQEGSLSARWPVAGQIRATQQQHDFTLTQTADGRSSFDWLTGSSTDPLVDHTSPSSDSLLFAHKRSERLQRAPLKSVGPDFGKKRLGLPGDEVDNKVKGAAGRTDLLRLRRRFMRDQEKLSLMYARKGVAEQKREKEIKSELKMKQDAQVVLYRSYRHGDLPDIQIKHSSLITPLQAVAQRDPIIAKQLFSSLFSGILKEMDKFKTLSEKNNITQKLLQDFNRFLNTTFSFFPPFVSCIQDISCQHAALLSLDPAAVSAGCLASLQQPVGIRLLEEALLRLLPAELPAKRVRGKARLPPDVLRWVELAKLYRSIGEYDVLRGIFTSEIGTKQITQSALLAEARSDYSEAAKQYDEALNKQDWVDGEPTEAEKDFWELASLDCYNHLAEWKSLEYCSTASIDSENPPDLNKIWSEPFYQETYLPYMIRSKLKLLLQGEADQSLLTFIDKAMHGELQKAILELHYSQELSLLYLLQDDVDRAKYYIQNGIQSFMQNYSSIDVLLHQSRLTKLQSVQALTEIQEFISFISKQGNLSSQVPLKRLLNTWTNRYPDAKMDPMNIWDDIITNRCFFLSKIEEKLTPLPEDNSMNVDQDGDPSDRMEVQEQEEDISSLIRSCKFSMKMKMIDSARKQNNFSLAMKLLKELHKESKTRDDWLVSWVQSYCRLSHCRSRSQGCSEQVLTVLKTVSLLDENNVSSYLSKNILAFRDQNILLGTTYRIIANALSSEPACLAEIEEDKARRILELSGSSSEDSEKVIAGLYQRAFQHLSEAVQAAEEEAQPPSWSCGPAAGVIDAYMTLADFCSFKDTSTGHKNKEFVARIKSKLDQGGVIQDFINALDQLSNPELLFKDWSNDVRAELAKTPVNKKNIEKMYERMYAALGDPKAPGLGAFRRKFIQTFGKEFDKHFGKGGSKLLRMKLSDFNDITNMLLLKMNKDSKPPGNLKECSPWMSDFKVEFLRNELEIPGQYDGRGKPLPEYHVRIAGFDERVTVMASLRRPKRIIIRGHDEREHPFLVKGGEDLRQDQRVEQLFQVMNGILAQDSACSQRALQLRTYSVVPMTSRLGLIEWLENTVTLKDLLLNTMSQEEKAAYLSDPRAPPCEYKDWLTKMSGKHDVGAYMLMYKGANRTETVTSFRKRESKVPADLLKRAFVRMSTSPEAFLALRSHFASSHALICISHWILGIGDRHLNNFMVAMETGGVIGIDFGHAFGSATQFLPVPELMPFRLTRQFINLMLPMKETGLMYSIMVHALRAFRSDPGLLTNTMDVFVKEPSFDWKNFEQKMLKKGGSWIQEINVAEKNWYPRQKICYAKRKLAGANPAVITCDELLLGHEKAPAFRDYVAVARGSKDHNIRAQEPESGLSEETQVKCLMDQATDPNILGRTWEGWEPWM"

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
# Each source of data has a location, such as a file on disk or a database
# entry, and a type. In this example we used EM density data, which we'll
# say lives in the EMDB database:
l = ihm.location.PDBLocation('5luq')
pdb_dataset = ihm.dataset.PDBDataset(l)

l_bsp = ihm.location.InputFileLocation("../data/BSP-350X-SEC.csv", details="BSP Crosslink File")
l_dss = ihm.location.InputFileLocation("../data/DSS-350X-SEC.csv", details="DSS Crosslink File")
l_dsg = ihm.location.InputFileLocation("../data/DSG-Crosslinks-2.csv", details="DSG Crosslink File")

bsp_dataset = ihm.dataset.CXMSDataset(l_bsp)
dss_dataset = ihm.dataset.CXMSDataset(l_dss)
dsg_dataset = ihm.dataset.CXMSDataset(l_dsg)


# Next, define the entities for each unique sequence in the system
entity = ihm.Entity(sequence, description='DNA-PKcs')
system.entities.append(entity)

# Next, we define asymmetric units for everything we modeled.
# These roughly correspond to chains in a traditional PDB file. Multiple
# asymmetric units may map to the same entity (for example if there are
# several copies of a given protein). Parts of the system that were seen in
# an experiment but were not modeled are represented as entities to which no
# asymmetric units map.
asym = ihm.AsymUnit(entity, details='Threading Model')
system.asym_units.append(asym)

# Next, we group asymmetric units (and/or entities) into assemblies.
# Here, we'll define an assembly of everything that we modeled, plus
# two subassemblies (of the subunits) that the SAXS data applies to:
assembly = ihm.Assembly([asym], name='Threading Model')


# Define how the system was represented. Multiple representations of the
# system are possible, and can overlap. Here we'll say we represent A
# atomically as a rigid body and B as 3 flexible coarse-grained spheres:
rep = ihm.representation.Representation(
        [ihm.representation.FeatureSegment(asym, rigid=True, primitive='sphere', count=4128, starting_model=pdb_dataset)])

# Set up restraints on the system. First, two on the subunits that use
# the SAXS data; we'll say we used the FoXS software to do this fit:

dsg = ihm.ChemDescriptor('DSG', chemical_name='disuccinimidyl glutarate',
              smiles='C1CC(=O)N(C1=O)OC(=O)CCCC(=O)ON2C(=O)CCC2=O',
              inchi='1S/C13H14N2O8/c16-8-4-5-9(17)14(8)22-12(20)2-1-3-13(21)23-15-10(18)6-7-11(15)19/h1-7H2',
              inchi_key='LNQHREYHFRFJAU-UHFFFAOYSA-N')

bsp = ihm.ChemDescriptor('BSP', chemical_name='Bis(succinimidyl) penta(ethylene glycol)',
              smiles='C1C(N(C(C1)=O)OC(CCOCCOCCOCCOCCOCCC(=O)ON2C(CCC2=O)=O)=O)=O',
              inchi='InChI=1S/C22H32N2O13/c25-17-1-2-18(26)23(17)36-21(29)5-7-31-9-11-33-13-15-35-16-14-34-12-10-32-8-6-22(30)37-24-19(27)3-4-20(24)28/h1-16H2',
              inchi_key='FTYUGLBWKRXQBD-UHFFFAOYSA-N')

dsg_restraint = ihm.restraint.CrossLinkRestraint(dsg_dataset, dsg)
dss_restraint = ihm.restraint.CrossLinkRestraint(dsg_dataset, ihm.cross_linkers.dss)
bsp_restraint = ihm.restraint.CrossLinkRestraint(dsg_dataset, bsp)
system.restraints.extend((dsg_restraint, dss_restraint, bsp_restraint))


# Now we add information about how the modeling was done by defining one
# or more protocols. Here we'll say we did simple Monte Carlo on the entire
# system using all of the experimental data:
all_datasets = ihm.dataset.DatasetGroup((pdb_dataset, dsg_dataset,
                                         dss_dataset, bsp_dataset))
protocol = ihm.protocol.Protocol(name='Modeling')
protocol.steps.append(ihm.protocol.Step(
                        assembly=assembly,
                        dataset_group=all_datasets,
                        method='Enumeration',
                        name='Production sampling',
                        num_models_begin=0,
                        num_models_end=1000, multi_scale=False))

# Finally we can add coordinates for the deposited models. Typically these
# will be stored in our own software's data structures somewhere (for this
# example in simple lists 'atoms' and 'spheres'):
atoms = [('A', 1, 'C', 'CA', 1., 2., 3.),
         ('A', 2, 'C', 'CA', 4., 5., 6.),
         ('A', 3, 'C', 'CA', 7., 8., 9.)]
spheres = [('B', 1, 2, 1., 2., 3., 1.2),
           ('B', 3, 4, 4., 5., 6., 1.2),
           ('B', 5, 6, 7., 8., 9., 1.2)]

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

    def get_atoms(self):
        # Use BioPython to read the structure from a PDB file, and then yield
        # a set of ihm.model.Atom objects
        p = Bio.PDB.PDBParser()
        s = p.get_structure('rep', self.file_name)
        for model in s:
            for nchain, chain in enumerate(model):
                asym = self.asym_units[nchain]
                for nres, residue in enumerate(chain):
                    for atom in residue:
                        coord = atom.get_vector()
                        yield ihm.model.Atom(asym_unit=asym, seq_id=nres+1,
                                atom_id=atom.get_id(), x=coord[0], y=coord[1],
                                z=coord[2], type_symbol=atom.element)



models = []
for pdb_file in glob.glob("./bsms_pdbs/*00.pdb"):
    print pdb_file
    m = Model(assembly = assembly, protocol=protocol, representation=rep,
        file_name=pdb_file, asym_units=[asym])
    s = ihm.model.State([ihm.model.ModelGroup([m])],
                        type='Threading Ensemble',
                        name='Threading Ensemble Solution')

    models.append(s)

mg = ihm.model.ModelGroup(models)
me = ihm.model.Ensemble(mg, len(models),
    post_process="filtering")

# Groups are then placed into states, which can in turn be grouped. In this
# case we have only a single state:
state = ihm.model.State([mg])
system.state_groups.append(ihm.model.StateGroup([state]))

# Once the system is complete, we can write it out to an mmCIF file:
with open('output.cif', 'w') as fh:
    ihm.dumper.write(fh, [system])