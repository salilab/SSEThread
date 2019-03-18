from __future__ import print_function
import IMP
import numpy
import IMP.pmi.tools
import IMP.atom
import IMP.core
import IMP.threading
import IMP.pmi.samplers
import json
import sys


# These arguments facilitate the multi-threading of enumeration runs
# They divide the total list of models into n_per_thread chunks and assigns
# chunk n_thread to this instance.
n_per_thread = int(sys.argv[1])
n_thread = int(sys.argv[2])
perm_number = int(sys.argv[3]) # Permutation of the SE sequence

# Hard-coded in permutations for a sequence of three
perms = [(0,1,2), (0,2,1), (1,0,2), (1,2,0), (2,1,0), (2,0,1)]
this_perm = perms[perm_number]

"""
 Normalised distances between loop lengths. The array index corresponds to loop length. The values are calculated for upto 30 residues.
The values are separately calculated for different secondary structural elements flanking the loop length. i.e there are 4 different values for
helix-loop-helix (nhh_mean), sheet-loop-helix(nsh_mean), helix-loop_sheet(nhs_mean) and sheet_loop_sheet(nss_mean).
Please note that the values are normalised against the loop_length
"""

nsh_mean = [0, 3.809, 3.137, 2.818, 2.482, 2.154, 1.928, 1.749, 1.67, 1.531, 1.428, 1.377, 1.282, 1.261, 1.203,
                1.135, 1.045, 1.004, 1.02, 0.977, 0.928, 0.865, 0.834, 0.811, 0.756, 0.761, 0.749, 0.777, 0.74, 0.655,
                0.648]
nhs_mean = [0, 3.809, 3.137, 2.818, 2.482, 2.154, 1.928, 1.749, 1.67, 1.531, 1.428, 1.377, 1.282, 1.261, 1.203,
                1.135, 1.045, 1.004, 1.02, 0.977, 0.928, 0.865, 0.834, 0.811, 0.756, 0.761, 0.749, 0.777, 0.74, 0.655,
                0.648]
nhh_mean = [0, 3.81, 3.036, 2.836, 2.511, 2.275, 2.178, 2.026, 1.876, 1.835, 1.669, 1.658, 1.666, 1.625, 1.53,
                1.445, 1.374, 1.292, 1.212, 1.164, 1.133, 1.049, 1.043, 1.074, 0.977, 0.965, 0.938, 0.868, 0.824, 0.805,
                0.788]
nss_mean = [0, 3.81, 3.19, 1.846, 1.607, 1.274, 1.14, 1.139, 1.198, 1.177, 1.115, 1.029, 1.048, 0.935, 0.91, 0.908,
                0.85, 0.83, 0.852, 0.849, 0.761, 0.722, 0.742, 0.684, 0.677, 0.611, 0.587, 0.596, 0.565, 0.576, 0.532]


# Standard deviations of the above normalised distances

hh_std = [0, 0.027, 0.284, 0.397, 0.441, 0.483, 0.499, 0.504, 0.537, 0.534, 0.538, 0.545, 0.507, 0.494, 0.468,
              0.447, 0.428, 0.439, 0.415, 0.432, 0.392, 0.382, 0.38, 0.401, 0.381, 0.38, 0.317, 0.328, 0.304, 0.318,
              0.273]
ss_std = [0, 0.027, 0.313, 0.293, 0.469, 0.419, 0.474, 0.49, 0.505, 0.447, 0.501, 0.475, 0.479, 0.417, 0.451, 0.416,
              0.373, 0.395, 0.47, 0.418, 0.36, 0.349, 0.359, 0.312, 0.302, 0.281, 0.279, 0.264, 0.259, 0.346, 0.257]
sh_std = [0, 0.067, 0.278, 0.361, 0.418, 0.45, 0.448, 0.455, 0.436, 0.452, 0.438, 0.416, 0.407, 0.402, 0.411, 0.405,
              0.381, 0.378, 0.373, 0.36, 0.372, 0.338, 0.322, 0.308, 0.285, 0.289, 0.296, 0.298, 0.294, 0.286, 0.208]
hs_std = [0, 0.067, 0.278, 0.361, 0.418, 0.45, 0.448, 0.455, 0.436, 0.452, 0.438, 0.416, 0.407, 0.402, 0.411, 0.405,
              0.381, 0.378, 0.373, 0.36, 0.372, 0.338, 0.322, 0.308, 0.285, 0.289, 0.296, 0.298, 0.294, 0.286, 0.208]

def randomize_ses_start_residues(ses, seq_length):
    # given a number of SEs, randomize the starting residues.
    import random
    avail_seq = range(1,seq_length+2)
    random.shuffle(ses)
    subseqs = [avail_seq]
    for s in ses:
        s_len = s.get_length()
        available_spots = []
        for ss in subseqs:
            if len(ss) > s_len:
                available_spots += ss[:-(s_len+2)]

        #print available_spots

        # pick an available spot
        spot = random.choice(available_spots)

        s.set_start_res_key(float(spot))

        for i in range(len(subseqs)):
            if spot in subseqs[i]:
                ss1 = range(subseqs[i][0],spot)
                ss2 = range(spot+s_len,subseqs[i][-1])

                prev_si = subseqs[i]
                del subseqs[i]
                if len(ss1) > 0: 
                    subseqs.append(ss1)
                if len(ss2) > 0:
                    subseqs.append(ss2)
                break


        #print spot, subseqs        

def read_psipred_ss2(ss2_file, ps):
    '''
    Reads a psipred .ss2 file and returns a list of 
    SecondaryStructureResidue particles that contains the 
    propensity values for Helix, Sheet and Loop that can be 
    inputted used for a IMP.threading.SecondaryStrctureRestraint object
    '''

    f = open(ss2_file, 'r')

    ss_probs = []
    for line in f.readlines():
        fields = line.strip().split()
        try: 
            int(fields[0])
            ss_probs.append((fields[4], fields[5], fields[3]))
        except:
            continue

    for i in range(len(ss_probs)):
        IMP.atom.SecondaryStructureResidue.setup_particle(ps[i], float(ss_probs[i][0]), float(ss_probs[i][1]), float(ss_probs[i][2]))

    return ps 

def sort_ses(ses):
    '''
    Given a list of structural elements, sort them by increasing residue number and return that list
    '''
    res = sorted([(s.get_first_residue_number(), s) for s in ses], key=lambda x: x[0])
    return [x[1] for x in res]

def setup_terminal_element(root_hier, coordinate, resid):
    '''
    Sets up a dummy StructureElement for a residue fixed in sequence space
    
    When building a threading model where some of the structure of one or more proteins is sequence assigned,
    we set up a StructureElement object for the last and first structured residues in the sequence loop(s) we are considering.
    These StructureElements contain the residue coordinate and number for use in computing the StructureElementConnectivityRestraint
    between this residue and the whatever sampled StructureElement is adjacent. This function also adds it to the system.

    @param root_hier : The root hierarchy of the system
    @param coordinate : The C_alpha coordinate of the endpoint residue
    @param resid : THe residue number of the endpoint residue
    
    Returns the generated StructureElement
    '''
    pi = IMP.Particle(root_hier.get_model())
    h = IMP.atom.Hierarchy.setup_particle(pi)
    np = IMP.Particle(root_hier.get_model())
    hp = IMP.atom.Hierarchy.setup_particle(np)
    xyz = IMP.core.XYZ.setup_particle(np)
    xyz.set_coordinates(coordinate)
    h.add_child(hp)

    # Set up Particle as StructureElement. Since it is only one residue, length=1, polarity does not matter and Offset=0.

    IMP.threading.StructureElement.setup_particle(root_hier.get_model(), pi.get_index(), resid, 1, 1, 0)
    
    # We do not optimize these keys because the residue assignment in sequence is known.
    se = IMP.threading.StructureElement(pi)
    se.set_keys_are_optimized(False)

    return se


def setup_structural_element(root_hier, element, max_translation=1):
    '''
    Set up a StructureElement for a given set of parameters, defined in element.
    element is a tuple (int, int, char) where the first int is the first residue in the PDB structure of this element
    and the second int is the number of residues in the element. char is either "H" for helix, "S" for sheet of "C" for coil.

    This function requires that a chain with ID "A" exists in root_hier that contains coordinates for these residue numbers.

    @param root_hier : The system root hierarchy
    @param element : The tuple contining the structure element information (start_res, length, SS assignment)
    @param max_translation : an integer that defines the maximum delta of the start_residue key per Monte Carlo move 
    '''

    # Add conditional for the existence of Chain A and coordinates for this set of residues

    # Grab particles corresponding to these residues from structure Chain A
    pis = IMP.atom.Selection(root_hier, chain_id="A", 
                residue_indexes=range(element[0],element[0]+element[1]),
                atom_type=IMP.atom.AT_CA).get_selected_particles()

    # Create new particle for this SE and set up as a Hierarchy
    pi = IMP.Particle(root_hier.get_model())
    h = IMP.atom.Hierarchy.setup_particle(pi)

    # Get XYZs from particles grabbed above and add to SE Hierarchy
    xyz = []
    i = 0
    for p in pis:
        np = IMP.Particle(root_hier.get_model())
        hp = IMP.atom.Hierarchy.setup_particle(np)
        xyz = IMP.core.XYZR.setup_particle(np)
        xyz.set_coordinates(IMP.core.XYZ(p).get_coordinates())
        xyz.set_radius(1.0)
        IMP.atom.Mass.setup_particle(np, 1.0)
        h.add_child(hp)
    
    # Set up StructureElement
    IMP.threading.StructureElement.setup_particle(root_hier.get_model(), pi.get_index(), element[0], 1, element[1], 0)

    # Set up the correct Secondary Structure assignment for this Element. Modeled as propensity = 1 for the SS assignment.

    if element[2]=="H":
        IMP.atom.SecondaryStructureResidue.setup_particle(pi, 1, 0, 0)
    elif element[2]=="S":        
        IMP.atom.SecondaryStructureResidue.setup_particle(pi, 0, 1, 0)
    elif element[2]=="C":
        IMP.atom.SecondaryStructureResidue.setup_particle(pi, 0, 0, 1)
    else:
        raise Exception("Secondary structure designation must be H, S or C.", element[2], "was provided")
    
    se = IMP.threading.StructureElement(pi)
    se.set_keys_are_optimized(True)

    return se

def setup_conditional_pair_restraint(p1, p2, length, xl_slope, constant):
    '''
    Sets up a ConditionalPairRestraint on two residues in sequence. For crosslinking or
    other pairwise distance restraints. Here, the restraint is set up as an upper harmonic 
    restraint with center at xl_length and stiffness xl_slope. 

    @param p1 : residue 1 particle
    @param p2 : residue 2 particle
    @param length : length in angstroms of the center of the upperharmonic (e.g. XL length plus side chains + tolerance)
    @param xl_slope : Spring constant of the upperharmonic function. (0.1 is a weak restraint. 1 is fairly strong)
    @param constant : The number of standard deviations above the mean to evaluate the loop length statistical potential. Can be negative if you want, but that would be a bit aggressive.
    '''
    r = IMP.threading.LoopPairDistanceRestraint(m, IMP.core.HarmonicUpperBound(length, xl_slope), 
        p1, p2, constant)

    return r

def setup_pair_restraint(p1, p2, length):
    r = IMP.core.DistanceRestraint(m, IMP.core.HarmonicUpperBound(length, xl_slope), 
        p1, p2)
    return r

def setup_length_restraint(s, length_slope):
    '''
    Setup a restraint that biases the StructureElement length key towards
    the length of the number of coordinates in the SE.  Implemented as a 
    Linear restraint on the difference between the two values.  
    
    @param s : The StructureElement object
    @param length_slope : The slope of the linear restraint. 
    '''
    uf = IMP.core.Linear(s.get_number_of_coordinates(), -1*length_slope)
    sf = IMP.core.AttributeSingletonScore(uf, IMP.FloatKey("length"))
    r = IMP.core.SingletonRestraint(m, sf, s.get_particle())
    return r

def add_SECR(p1, p2, slope=1, n_sds=1.0):
    '''
    Define a StructureElementConnectivityRestraint between two StructureElements

    @param p1 : ParticleIndex of first StructureElement (N-terminal)
    @param p2 : ParticleIndex of second StructureElement (C-terminal)
    @param slope : Spring constant of the harmonic upper bound
    @param n_sds : Number of standard deviations above the mean at which to compute the statistical potential distance
    '''
    r = IMP.threading.StructureElementConnectivityRestraint(m, IMP.core.HarmonicUpperBound(0, slope), p1, p2, n_sds, "")
    return r

def add_all_SECR(se_list, slope=1, dpr=1.0):
    '''
    Simple function to create many StructureElementConnectivityRestraints from a list of ParticleIndexes
    of ordered StructureElements
    '''
    print("Adding", len(se_list)-1, "StructureElementConnectivityRestraints.")
    SECR_restraints = []
    for i in range(len(se_list-1)):
        p1 = se_list[i].get_particle_index()
        p2 = se_list[i+1].get_particle_index()

        r = IMP.threading.StructureElementConnectivityRestraint(m, IMP.core.HarmonicUpperBound(0, slope), p1, p2, dpr, "")
        SECR_restraints.append(r)

    return SECR_restraints

def modify_all_SECR(se_list, rst_list):

    for i in range(len(rst_list)):
        p1 = se_list[i].get_particle_index()
        p2 = se_list[i+1].get_particle_index()
        rst_list.assign_particles(p1, p2)

    return SECR_restraints

#  Sequence of DNAPKcs
seq = "MAGSGAGVRCSLLRLQETLSAADRCGAALAGHQLIRGLGQECVLSSSPAVLALQTSLVFSRDFGLLVFVRKSLNSIEFRECREEILKFLCIFLEKMGQKIAPYSVEIKNTCTSVYTKDRAAKCKIPALDLLIKLLQTFRSSRLMDEFKIGELFSKFYGELALKKKIPDTVLEKVYELLGLLGEVHPSEMINNAENLFRAFLGELKTQMTSAVREPKLPVLAGCLKGLSSLLCNFTKSMEEDPQTSREIFNFVLKAIRPQIDLKRYAVPSAGLRLFALHASQFSTCLLDNYVSLFEVLLKWCAHTNVELKKAALSALESFLKQVSNMVAKNAEMHKNKLQYFMEQFYGIIRNVDSNNKELSIAIRGYGLFAGPCKVINAKDVDFMYVELIQRCKQMFLTQTDTGDDRVYQMPSFLQSVASVLLYLDTVPEVYTPVLEHLVVMQIDSFPQYSPKMQLVCCRAIVKVFLALAAKGPVLRNCISTVVHQGLIRICSKPVVLPKGPESESEDHRASGEVRTGKWKVPTYKDYVDLFRHLLSSDQMMDSILADEAFFSVNSSSESLNHLLYDEFVKSVLKIVEKLDLTLEIQTVGEQENGDEAPGVWMIPTSDPAANLHPAKPKDFSAFINLVEFCREILPEKQAEFFEPWVYSFSYELILQSTRLPLISGFYKLLSITVRNAKKIKYFEGVSPKSLKHSPEDPEKYSCFALFVKFGKEVAVKMKQYKDELLASCLTFLLSLPHNIIELDVRAYVPALQMAFKLGLSYTPLAEVGLNALEEWSIYIDRHVMQPYYKDILPCLDGYLKTSALSDETKNNWEVSALSRAAQKGFNKVVLKHLKKTKNLSSNEAISLEEIRIRVVQMLGSLGGQINKNLLTVTSSDEMMKSYVAWDREKRLSFAVPFREMKPVIFLDVFLPRVTELALTASDRQTKVAACELLHSMVMFMLGKATQMPEGGQGAPPMYQLYKRTFPVLLRLACDVDQVTRQLYEPLVMQLIHWFTNNKKFESQDTVALLEAILDGIVDPVDSTLRDFCGRCIREFLKWSIKQITPQQQEKSPVNTKSLFKRLYSLALHPNAFKRLGASLAFNNIYREFREEESLVEQFVFEALVIYMESLALAHADEKSLGTIQQCCDAIDHLCRIIEKKHVSLNKAKKRRLPRGFPPSASLCLLDLVKWLLAHCGRPQTECRHKSIELFYKFVPLLPGNRSPNLWLKDVLKEEGVSFLINTFEGGGCGQPSGILAQPTLLYLRGPFSLQATLCWLDLLLAALECYNTFIGERTVGALQVLGTEAQSSLLKAVAFFLESIAMHDIIAAEKCFGTGAAGNRTSPQEGERYNYSKCTVVVRIMEFTTTLLNTSPEGWKLLKKDLCNTHLMRVLVQTLCEPASIGFNIGDVQVMAHLPDVCVNLMKALKMSPYKDILETHLREKITAQSIEELCAVNLYGPDAQVDRSRLAAVVSACKQLHRAGLLHNILPSQSTDLHHSVGTELLSLVYKGIAPGDERQCLPSLDLSCKQLASGLLELAFAFGGLCERLVSLLLNPAVLSTASLGSSQGSVIHFSHGEYFYSLFSETINTELLKNLDLAVLELMQSSVDNTKMVSAVLNGMLDQSFRERANQKHQGLKLATTILQHWKKCDSWWAKDSPLETKMAVLALLAKILQIDSSVSFNTSHGSFPEVFTTYISLLADTKLDLHLKGQAVTLLPFFTSLTGGSLEELRRVLEQLIVAHFPMQSREFPPGTPRFNNYVDCMKKFLDALELSQSPMLLELMTEVLCREQQHVMEELFQSSFRRIARRGSCVTQVGLLESVYEMFRKDDPRLSFTRQSFVDRSLLTLLWHCSLDALREFFSTIVVDAIDVLKSRFTKLNESTFDTQITKKMGYYKILDVMYSRLPKDDVHAKESKINQVFHGSCITEGNELTKTLIKLCYDAFTENMAGENQLLERRRLYHCAAYNCAISVICCVFNELKFYQGFLFSEKPEKNLLIFENLIDLKRRYNFPVEVEVPMERKKKYIEIRKEAREAANGDSDGPSYMSSLSYLADSTLSEEMSQFDFSTGVQSYSYSSQDPRPATGRFRRREQRDPTVHDDVLELEMDELNRHECMAPLTALVKHMHRSLGPPQGEEDSVPRDLPSWMKFLHGKLGNPIVPLNIRLFLAKLVINTEEVFRPYAKHWLSPLLQLAASENNGGEGIHYMVVEIVATILSWTGLATPTGVPKDEVLANRLLNFLMKHVFHPKRAVFRHNLEIIKTLVECWKDCLSIPYRLIFEKFSGKDPNSKDNSVGIQLLGIVMANDLPPYDPQCGIQSSEYFQALVNNMSFVRYKEVYAAAAEVLGLILRYVMERKNILEESLCELVAKQLKQHQNTMEDKFIVCLNKVTKSFPPLADRFMNAVFFLLPKFHGVLKTLCLEVVLCRVEGMTELYFQLKSKDFVQVMRHRDDERQKVCLDIIYKMMPKLKPVELRELLNPVVEFVSHPSTTCREQMYNILMWIHDNYRDPESETDNDSQEIFKLAKDVLIQGLIDENPGLQLIIRNFWSHETRLPSNTLDRLLALNSLYSPKIEVHFLSLATNFLLEMTSMSPDYPNPMFEHPLSECEFQEYTIDSDWRFRSTVLTPMFVETQASQGTLQTRTQEGSLSARWPVAGQIRATQQQHDFTLTQTADGRSSFDWLTGSSTDPLVDHTSPSSDSLLFAHKRSERLQRAPLKSVGPDFGKKRLGLPGDEVDNKVKGAAGRTDLLRLRRRFMRDQEKLSLMYARKGVAEQKREKEIKSELKMKQDAQVVLYRSYRHGDLPDIQIKHSSLITPLQAVAQRDPIIAKQLFSSLFSGILKEMDKFKTLSEKNNITQKLLQDFNRFLNTTFSFFPPFVSCIQDISCQHAALLSLDPAAVSAGCLASLQQPVGIRLLEEALLRLLPAELPAKRVRGKARLPPDVLRWVELAKLYRSIGEYDVLRGIFTSEIGTKQITQSALLAEARSDYSEAAKQYDEALNKQDWVDGEPTEAEKDFWELASLDCYNHLAEWKSLEYCSTASIDSENPPDLNKIWSEPFYQETYLPYMIRSKLKLLLQGEADQSLLTFIDKAMHGELQKAILELHYSQELSLLYLLQDDVDRAKYYIQNGIQSFMQNYSSIDVLLHQSRLTKLQSVQALTEIQEFISFISKQGNLSSQVPLKRLLNTWTNRYPDAKMDPMNIWDDIITNRCFFLSKIEEKLTPLPEDNSMNVDQDGDPSDRMEVQEQEEDISSLIRSCKFSMKMKMIDSARKQNNFSLAMKLLKELHKESKTRDDWLVSWVQSYCRLSHCRSRSQGCSEQVLTVLKTVSLLDENNVSSYLSKNILAFRDQNILLGTTYRIIANALSSEPACLAEIEEDKARRILELSGSSSEDSEKVIAGLYQRAFQHLSEAVQAAEEEAQPPSWSCGPAAGVIDAYMTLADFCDQQLRKEEENASVIDSAELQAYPALVVEKMLKALKLNSNEARLKFPRLLQIIERYPEETLSLMTKEISSVPCWQFISWISHMVALLDKDQAVAVQHSVEEITDNYPQAIVYPFIISSESYSFKDTSTGHKNKEFVARIKSKLDQGGVIQDFINALDQLSNPELLFKDWSNDVRAELAKTPVNKKNIEKMYERMYAALGDPKAPGLGAFRRKFIQTFGKEFDKHFGKGGSKLLRMKLSDFNDITNMLLLKMNKDSKPPGNLKECSPWMSDFKVEFLRNELEIPGQYDGRGKPLPEYHVRIAGFDERVTVMASLRRPKRIIIRGHDEREHPFLVKGGEDLRQDQRVEQLFQVMNGILAQDSACSQRALQLRTYSVVPMTSRLGLIEWLENTVTLKDLLLNTMSQEEKAAYLSDPRAPPCEYKDWLTKMSGKHDVGAYMLMYKGANRTETVTSFRKRESKVPADLLKRAFVRMSTSPEAFLALRSHFASSHALICISHWILGIGDRHLNNFMVAMETGGVIGIDFGHAFGSATQFLPVPELMPFRLTRQFINLMLPMKETGLMYSIMVHALRAFRSDPGLLTNTMDVFVKEPSFDWKNFEQKMLKKGGSWIQEINVAEKNWYPRQKICYAKRKLAGANPAVITCDELLLGHEKAPAFRDYVAVARGSKDHNIRAQEPESGLSEETQVKCLMDQATDPNILGRTWEGWEPWM"


# Restraint Weights
#---------------------
length_slope = 0.1
semet_slope = 0.1
xl_constant = 1.0    # Number of SDs higher than mean to evaluate model XL distance
xl_slope = 0.05      # Stiffness of XL restraint upper-harmonic function
psipred_slope = 0.1

# Crosslinked residues manually enetered from crosslinking data
# Crosslink evaluation distances (third entry in each tuple) are XL length plus 20 angstroms
#-------------------------
xl_22 = [(2746,99,42),(2764,810,42),(2754,810,42),(2764,838,42),(2738,2003,42),(2694,2003,42),(2717,2107,42),(2738,2313,42),(2694,2445,42)]
xl_12 = [(2746,99,32),(2746,357,32),(2764,810,32),(2738,2227,32),(2746,2313,32)]
xl_8 = [(2746,357,27),(2738,2001,27),(2717,2107,27),(2738,2227,27),(2746,2313,27)]

# Data files
#--------------------------
pdbfile = "../data/5luq_A_CA.pdb"
psipred_file = "../data/psipred_dnapkcs.txt"

# Define Structural Elements
# These are the residues in the PDB that correspond to structural elements.
# (start_res, length, SSID)
#-------------------------
elements=[(2602,14,'H'),(2623,25,'H'), (2655,10,'H')]


# Define terminal residue coordinates and residue number
#----------------------------
terminal_residues = [[(27.365, -37.659, 9.827),2575],[(49.775, -42.147, 7.299),2774]]

# Residue range of the disordered area we wish to sample 
# --------------------------
disordered_range=(2577,2772)

#####################
######################################
# Here is where the work begins
######################################
#####################


##########
# Set Up System Representation
##########

m = IMP.Model()

print("Building Structure Chain and StructureElements")

###
# Set up Structure chain A from PDB file
# ------------------

root_hier = IMP.atom.read_pdb(pdbfile, m)

# Create StructureElement particles
se = []
for e in elements:
    se.append(setup_structural_element(root_hier, e))

# Set up terminal residues StructureElements
term_ses = []
for t in terminal_residues:
    term_ses.append(setup_terminal_element(root_hier, t[0], t[1]))


###
# Set up Sequence chain S from seq
# ---------------

print("Building Sequence Chain for", len(seq), "residues")

seq_chain = IMP.atom.Chain.setup_particle(IMP.Particle(m), "S")
root_hier.add_child(seq_chain)

res_particles = []

# Add a residue for all residues in the input sequence
for i in range(len(seq)):

    pr = IMP.Particle(m)
    res = IMP.atom.Residue.setup_particle(pr,
                                        IMP.pmi.tools.get_residue_type_from_one_letter_code(seq[i]),
                                        i+1)
    res_particles.append(res.get_particle())
    IMP.atom.Mass.setup_particle(res.get_particle(), IMP.atom.get_mass(res.get_residue_type()))
    IMP.core.XYZR.setup_particle(res.get_particle())

    # Check to see if this residue exists in Chain A
    sel = IMP.atom.Selection(root_hier.get_children()[0], residue_index=i).get_selected_particles()
    
    # If so, assign the coordinate to this sequence residue in chain S
    if len(sel) == 1:
        IMP.core.XYZ(pr).set_coordinates(IMP.core.XYZ(sel[0]).get_coordinates())
        IMP.core.XYZ(pr).set_coordinates_are_optimized(True)
    
    seq_chain.add_child(res)


#######################
# Set up Scoring Function
#######################

# List to hold all restraints
rests = []

# XL Restraint
# ---------
# Restraint is calculated using Crosslinking distance and LoopPairDistanceRestraint, which takes into account
# intervening residues if one or more crosslinking sites are in a loop with no coordiantes
# Create list of crosslinked sites only between residues that exist in the model

xl_sites = []
for xl in xl_22 + xl_12 + xl_8:
    p = IMP.atom.Selection(root_hier, residue_index=xl[1]).get_selected_particles()
    if len(p) > 0:
        xl_sites.append(xl)

    # TODO: Add warning/exception for too many/no particles in selection

xlr = []

print("Adding", len(xl_sites), "crosslinking restraints.")
for xl in xl_sites:
    p1 = IMP.atom.Selection(root_hier, chain_id='S', residue_index=xl[0]).get_selected_particle_indexes()[0]
    p2 = IMP.atom.Selection(root_hier, chain_id='S', residue_index=xl[1]).get_selected_particle_indexes()[0]

    r = setup_conditional_pair_restraint(p1, p2, xl[2], xl_slope, xl_constant)
    pis = r.get_sequence_residue_particles()

    r1 = IMP.atom.Residue(m, pis[0])
    r2 = IMP.atom.Residue(m, pis[1])
    print("XL Setup between residues: ", xl[0], xl[1], " at a length of ", xl[2], "|", r1, r2)
    #print r.evaluate(False)
    rests.append(r)
    xlr.append(r)


# Completeness Restraint
# ---------
# Restraint on the number of StructureElement coordinates modeled

for s in se:
    r = setup_length_restraint(s, length_slope=length_slope)
    rests.append(r)


# Loop Length distances
# ---------
# Structural domain : sequence of SEs.  Use the sequence to 
# add the StructureElementConnectivity restraint.  model distance = last coord of SE1 - first coord of SE2.
# evaluated distance = (first resid of SE1 - last resid of SE2) * res_dist

# TODO: Automate this into a single-liner

# Define StructureElement pairs to add to connectivity restraint
# Assign their start_residues as per the permutation of this thread
se_pairs = []
se_pairs.append((term_ses[0].get_particle_index(), se[this_perm[0]].get_particle_index()))
se_pairs.append((se[this_perm[0]].get_particle_index(), se[this_perm[1]].get_particle_index()))
se_pairs.append((se[this_perm[1]].get_particle_index(), se[this_perm[2]].get_particle_index()))
se_pairs.append((se[this_perm[2]].get_particle_index(), term_ses[-1].get_particle_index()))


secrs = []
print ("Adding", len(se_pairs), "StructureElementConnectivityRestraints.")
for s in se_pairs[::-1]:
    r = add_SECR(s[0], s[1])
    secrs.append(r)
    rests.append(r)

# PSIPRED Restraint
# ---------
# A restraint on the overlap between the DSSP defined secondary structure of the StructureElement
# and the PSIPRED propensities of that residue in sequence

read_psipred_ss2(psipred_file, res_particles)

print("Psipred information added to Chain S")
print("SecondaryStructureParsimonyRestraint added")

se_part_indexes = [s.get_particle_index() for s in se]
r = IMP.threading.SecondaryStructureParsimonyRestraint(m, se_part_indexes, seq_chain.get_particle_index(), 1.0)
rests.append(r)

# Create a scoring function for all restraints
#-----------------------------------
sf = IMP.core.RestraintsScoringFunction(rests)


#########################
# Sampling Setup
#########################

print("Building Monte Carlo sampling objects")

# Define MonteCarlo sampling object and add scoring function
#--------------------
mc = IMP.core.MonteCarlo(m)
mc.set_scoring_function(sf)
mc.set_return_best(False)


# Set up StructureElement Samplers
#-------------------

sems = []

for s in se:
    coordinates = s.get_coordinates()
    sem = IMP.threading.StructureElementMover(m, s.get_particle_index(), seq_chain.get_particle())
    
    # Copy coordinates from StructureElement to Sequence chain
    sem.transform_coordinates()

    sems.append(sem)

    # Add the mover for sampling
    mc.add_mover(sem)


def coolate_output(se, sf):
    outlist = [[s.get_all_key_values() for s in se], sf.evaluate(False)]
    return outlist

def run_one_sim(mc, n_equil_frames, n_prod_frames):
    mc.set_return_best(False)
    for f in range(n_equil_frames):
        mc.optimize(1)
    mc.set_return_best(True)
    for f in range(n_prod_frames):
        mc.optimize(10)

print("Total Score on initialization:", sf.evaluate(False), [(r.get_name(), r.evaluate(False)) for r in rests])

all_output=[]

def enumerate_start_resis_3(seq_bounds, se_lengths):
    # Given the lengths of three SEs and the residue bounds, enumerate the potential start residues sets for each
    # without overlaps.  The order of the SEs is not changed.

    # TODO: allow for disjoint residue ranges
    # TODO: Make general for N StructureElements

    import itertools

    seq_len = seq_bounds[1] - seq_bounds[0]

    # Get all permutations
    all_se_perms = list(itertools.permutations(se_lengths, len(se_lengths)))

    start_resis = []
    for i in range(seq_bounds[0], seq_bounds[1] + 1 -se_lengths[0]-se_lengths[1]-se_lengths[2]-6):
        for j in range(se_lengths[0]+i+2, seq_bounds[1]+1-se_lengths[1]-se_lengths[2]-4):
            for k in range(se_lengths[1]+j+2, seq_bounds[1]+1-se_lengths[2]-2):
                start_resis.append((i,j,k))

    return start_resis


def enumerate_start_resis_3_all(seq_bounds, se_lengths, force=False):
    # Given three SE's, enumerate the potential start residues.
    # Permute the order of the SEs 

    print("Setting up enumeration for", len(se_lengths), "SEs", "from residues", seq_bounds[0], "to", seq_bounds[1])
    
    if len(se_lengths)!=3:
        raise Exception("Must have only three SEs to enumerate with enumerate_start_residue_3_all().", len(se_lengths), "were passed")

    # First set is in the order given
    ses = enumerate_start_resis_3(seq_bounds, se_lengths)

    # TODO: Put these in an automatically generated loop

    # 0 2 1
    new_se_lengths = [se_lengths[0], se_lengths[2], se_lengths[1]]
    new_ses = enumerate_start_resis_3(seq_bounds, new_se_lengths)
    # Now swap 1 and 2 in the output
    for i in new_ses:
        ses.append((i[0], i[2], i[1]))

    # 1 0 2
    new_se_lengths = [se_lengths[1], se_lengths[0], se_lengths[2]]
    new_ses = enumerate_start_resis_3(seq_bounds, new_se_lengths)
    # Now swap 1 and 2 in the output
    for i in new_ses:
        ses.append((i[1], i[0], i[2]))

    # 1 2 0
    new_se_lengths = [se_lengths[1], se_lengths[2], se_lengths[0]]
    new_ses = enumerate_start_resis_3(seq_bounds, new_se_lengths)
    # Now swap 1 and 2 in the output
    for i in new_ses:
        ses.append((i[1], i[2], i[0]))

    # 2 1 0
    new_se_lengths = [se_lengths[2], se_lengths[1], se_lengths[0]]
    new_ses = enumerate_start_resis_3(seq_bounds, new_se_lengths)
    # Now swap 1 and 2 in the output
    for i in new_ses:
        ses.append((i[2], i[1], i[0]))

    # 2 0 1
    new_se_lengths = [se_lengths[2], se_lengths[0], se_lengths[1]]
    new_ses = enumerate_start_resis_3(seq_bounds, new_se_lengths)
    # Now swap 1 and 2 in the output
    for i in new_ses:
        ses.append((i[2], i[0], i[1]))

    return ses


def enumerate_start_resis_3_eval(start_resis, ses, sems, sf, f):
    '''
    Given a set of models (start_resis), evaluate each of these models
    on the set of StructureElements (ses) with associated StructureElementMovers (sems)
    using scoring function (sf) 
    '''
    
    out_dict={}
    out_dict['models'] = []

    # TODO: Change to simply adding Chain A and grab the SEs and SEMs from that rather than passing them

    if len(ses) != len(start_resis[0]):
        raise Exception("Different number of StructureElements,", len(ses),"and model values,", len(start_resis[0]))
    if len(ses) != len(sems):
        raise Exception("Different number of StructureElements,", len(ses),"and StructureElementMovers,", len(sems))

    for i in range(len(start_resis)):
        # Apply start residues to each SE
        for n in range(len(start_resis[i])):
            sr = start_resis[i][n]
            ste = ses[n]
            sem = sems[n]
            ste.set_start_res_key(sr)
            sem.transform_coordinates()
        
        # log the model, frame and value of all scoring functions 
        new_mod = {}
        new_mod['frame'] = i
        new_mod['model'] = [s.get_all_key_values() for s in se]
        new_mod['restraints'] = {}
        tscore = 0
        for r in rests:
            score = r.evaluate(False)
            tscore+=score
            new_mod['restraints'][r.get_name()] = r.evaluate(False)
        new_mod['score'] = tscore

        out_dict['models'].append(new_mod) 

        # Print every 100th step just to know that we're actually moving forward
        if i%100==0: 
            print(str(i)+"th step", tscore) #sf.evaluate(False), start_resis[i]
    return out_dict 


# Define Sampling Space
# ------------------
# Enumerate the possible start_residue sets within the disordered_range

se_lengths = [e[1] for e in elements]

print("SE:", se_lengths)

# Re-order list of SE lengths
our_se_lengths = [se_lengths[this_perm[0]], se_lengths[this_perm[1]], se_lengths[this_perm[2]]]

# Enumerate possible models for three SEs
srs = enumerate_start_resis_3_all((disordered_range[0], disordered_range[1]-our_se_lengths[-1]), our_se_lengths)

print("Total models to evaluate:", len(srs))

# Indexes to grab the set of models evaluated in this thread
first = n_thread * n_per_thread
last = first + n_per_thread
srs_for_us = srs[first:last]

new_srs = []
for s in srs_for_us:
   new_srs.append((s[this_perm.index(0)], s[this_perm.index(1)],s[this_perm.index(2)]))

print("Thread", n_thread, "is evaluating models", first, "to", last)


#######################
# Define Output and begin evaluating models
#######################

import time
a = time.time()
f = open("enumerate_multi_"+str(n_thread)+"_"+str(n_per_thread)+".dat", "w")

# This function evaluates the set of models
out_dict = enumerate_start_resis_3_eval(srs_for_us, se, sems, sf, f)
a1 = time.time()
print("Final "+str(len(out_dict['models'])) +" models added")
print("Total time ", a1-a, "Time per model =", (a1-a)/len(srs_for_us))

# Write to output file
json.dump(out_dict, f) 
out_dict={}
IMP.atom.destroy(root_hier)

# TODO: Fix bug where entire IMP hierarchy is not destroyed


