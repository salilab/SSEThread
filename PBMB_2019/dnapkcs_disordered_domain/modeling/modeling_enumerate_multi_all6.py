import IMP
import numpy
import IMP.pmi.tools
import IMP.atom
import IMP.core
import IMP.threading
import IMP.pmi.samplers
import json
import time
import sys

# NMR structure of the FF domain L24A mutant's folding transition state
print sys.argv
n_per_thread = int(sys.argv[1])
n_thread = int(sys.argv[2])
perm_number = int(sys.argv[3])

perms = [(0,1,2), (0,2,1), (1,0,2), (1,2,0), (2,1,0), (2,0,1)]

this_perm = perms[perm_number]

print n_per_thread, n_thread, n_per_thread * n_thread, perm_number

init=time.time()
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



def enumerate_start_resis_3(seq_bounds, se_lengths, force=False):
    # Given three SE's, enumerate the potential start residues
    import itertools

    if se_lengths > 5 and force is True:
        print "Are you sure you want to enumerate this many SEs?"
        print "If so, set force=True"
        exit()

    seq_len = seq_bounds[1] - seq_bounds[0]

    print "Setting up enumeration for", len(se_lengths), "SEs", "from residues", seq_bounds[0], "to", seq_bounds[1]

    all_se_perms = list(itertools.permutations(se_lengths, len(se_lengths)))

    start_resis = []

    for i in range(seq_bounds[0], seq_bounds[1] + 1 -se_lengths[0]-se_lengths[1]-se_lengths[2]-6):
        for j in range(se_lengths[0]+i+2, seq_bounds[1]+1-se_lengths[1]-se_lengths[2]-4):
            for k in range(se_lengths[1]+j+2, seq_bounds[1]+1-se_lengths[2]-2):
                start_resis.append((i,j,k))


    return start_resis

def enumerate_start_resis_3_eval(start_resis, ses, sems, sf, f):

    out_dict={}
    out_dict['models'] = []

    print "Let's do", len(start_resis), "calculations!"

    for i in range(len(start_resis)):
        for n in range(len(start_resis[i])):
            sr = start_resis[i][n]
            ste = ses[n]
            sem = sems[n]
            ste.set_start_res_key(sr)
            sem.transform_coordinates()
 
        new_mod = {}
        new_mod['frame'] = i
        new_mod['model'] = [s.get_all_key_values() for s in ses]
        new_mod['restraints'] = {}
        tscore = 0
        #print "MOD: ", new_mod['model']
        for r in rests:
            score = r.evaluate(False)
            tscore+=score

            new_mod['restraints'][r.get_name()] = r.evaluate(False)
        new_mod['score'] = tscore

        out_dict['models'].append(new_mod) 
        if i%100==0: 
            print str(i)+"th step", tscore #sf.evaluate(False), start_resis[i]
      
    return out_dict 


def randomize_ses_start_residues(ses, seq_length):
    # given a number of SEs, randomize the start residues.
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

    f = open(ss2_file, 'r')

    ss_probs = []
    for line in f.readlines():
        fields = line.strip().split()
        try: 
            int(fields[0])
            ss_probs.append((fields[4], fields[5], fields[3]))
        except:
            continue

    #if len(ss_probs) != len(ps):
    #    raise Exception("SS file does not have same # of residues as protein")

    for i in range(len(ss_probs)):
        #print ps[i], ss_probs[i][0], ss_probs[i][1], ss_probs[i][2], ss_probs[i][0]
        IMP.atom.SecondaryStructureResidue.setup_particle(ps[i], float(ss_probs[i][0]), float(ss_probs[i][1]), float(ss_probs[i][2]))

    return ps 

def sort_ses(ses):
    # Given a list of structural elements, sort them by increasing first residue
    res = sorted([(s.get_first_residue_number(), s) for s in ses], key=lambda x: x[0])
    #print res
    return [x[1] for x in res]

def setup_terminal_element(root_hier, coordinate, resid):
    pi = IMP.Particle(root_hier.get_model())
    h = IMP.atom.Hierarchy.setup_particle(pi)
    np = IMP.Particle(root_hier.get_model())
    hp = IMP.atom.Hierarchy.setup_particle(np)
    #rp = IMP.atom.Residue.setup_particle(np)
    xyz = IMP.core.XYZ.setup_particle(np)
        #m = IMP.atom.Mass.setup_particle(root_hier.get_model(), np)
    xyz.set_coordinates(coordinate)
    h.add_child(hp)

    IMP.threading.StructureElement.setup_particle(root_hier.get_model(), pi.get_index(), resid, 1, 1, 0)


    se = IMP.threading.StructureElement(pi)
    se.set_keys_are_optimized(False)

    return se


def setup_structural_element(root_hier, element, max_translation=1):
    ## SETUP STRUCTURAL ELEMENTS
    # Get Particle Indexes for CA atoms 
    # from residue range element[0]:element[1]
    # From chain X
    # @param element : a tuple of form (start_res, length, 'H'/'E'/'C')
    #print(range(element[0],element[0]+element[1]))

    pis = IMP.atom.Selection(root_hier, chain_id="A", 
                residue_indexes=range(element[0],element[0]+element[1]),
                atom_type=IMP.atom.AT_CA).get_selected_particles()

    pi = IMP.Particle(root_hier.get_model())

    h = IMP.atom.Hierarchy.setup_particle(pi)

    # Get XYZs
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

    IMP.threading.StructureElement.setup_particle(root_hier.get_model(), pi.get_index(), element[0], 1, element[1], 0)

    # Set up this element as a helix
    # particle, Helix, sheet, coil
    IMP.atom.SecondaryStructureResidue.setup_particle(pi, 1, 0, 0)

    se = IMP.threading.StructureElement(pi)
    se.set_keys_are_optimized(True)
    se.set_max_res(200)

    return se

def setup_conditional_pair_restraint(p1, p2, length, xl_slope, constant):
    r = IMP.threading.LoopPairDistanceRestraint(m, IMP.core.HarmonicUpperBound(length, xl_slope), 
        p1, p2, constant)

    return r

def setup_pair_restraint(p1, p2, length):
    #dps = IMP.core.DistancePairScore(IMP.core.HarmonicUpperBound(length, xl_slope))
    r = IMP.core.DistanceRestraint(m, IMP.core.HarmonicUpperBound(length, xl_slope), 
        p1, p2)
    return r

def setup_length_restraint(s, length_slope):
    # Setup a restraint that biases the structural element towards
    # the length of the number of coordinates.

    # score = -log(#)
    uf = IMP.core.Linear(s.get_number_of_coordinates(), -1*length_slope)
    sf = IMP.core.AttributeSingletonScore(uf, IMP.FloatKey("length"))
    r = IMP.core.SingletonRestraint(m, sf, s.get_particle())
    print "SSR", r
    return r

def add_SECR(p1, p2, slope=1, dpr=1.0):
    r = IMP.threading.StructureElementConnectivityRestraint(m, IMP.core.HarmonicUpperBound(0, slope), p1, p2, dpr, "")
    return r

def add_all_SECR(se_list, slope=1, dpr=1.0):
    SECR_restraints = []
    for i in range(len(se_list-1)):
        p1 = se_list[i].get_particle_index()
        p2 = se_list[i+1].get_particle_index()

        r = IMP.threading.StructureElementConnectivityRestraint(m, IMP.core.HarmonicUpperBound(0, slope), p1, p2, dpr, "")
        SECR_restraints.append(r)

        print ">>>>>", p1, p2
        print r.get_number_of_residues(), r.get_max_distance, r.get_model_distance()
        print r.unprotected_evaluate(None)


    return SECR_restraints

def modify_all_SECR(se_list, rst_list):

    for i in range(len(rst_list)):
        p1 = se_list[i].get_particle_index()
        p2 = se_list[i+1].get_particle_index()
        rst_list.assign_particles(p1, p2)

    return SECR_restraints

#  Sequence of DNAPKcs
seq = "MAGSGAGVRCSLLRLQETLSAADRCGAALAGHQLIRGLGQECVLSSSPAVLALQTSLVFSRDFGLLVFVRKSLNSIEFRECREEILKFLCIFLEKMGQKIAPYSVEIKNTCTSVYTKDRAAKCKIPALDLLIKLLQTFRSSRLMDEFKIGELFSKFYGELALKKKIPDTVLEKVYELLGLLGEVHPSEMINNAENLFRAFLGELKTQMTSAVREPKLPVLAGCLKGLSSLLCNFTKSMEEDPQTSREIFNFVLKAIRPQIDLKRYAVPSAGLRLFALHASQFSTCLLDNYVSLFEVLLKWCAHTNVELKKAALSALESFLKQVSNMVAKNAEMHKNKLQYFMEQFYGIIRNVDSNNKELSIAIRGYGLFAGPCKVINAKDVDFMYVELIQRCKQMFLTQTDTGDDRVYQMPSFLQSVASVLLYLDTVPEVYTPVLEHLVVMQIDSFPQYSPKMQLVCCRAIVKVFLALAAKGPVLRNCISTVVHQGLIRICSKPVVLPKGPESESEDHRASGEVRTGKWKVPTYKDYVDLFRHLLSSDQMMDSILADEAFFSVNSSSESLNHLLYDEFVKSVLKIVEKLDLTLEIQTVGEQENGDEAPGVWMIPTSDPAANLHPAKPKDFSAFINLVEFCREILPEKQAEFFEPWVYSFSYELILQSTRLPLISGFYKLLSITVRNAKKIKYFEGVSPKSLKHSPEDPEKYSCFALFVKFGKEVAVKMKQYKDELLASCLTFLLSLPHNIIELDVRAYVPALQMAFKLGLSYTPLAEVGLNALEEWSIYIDRHVMQPYYKDILPCLDGYLKTSALSDETKNNWEVSALSRAAQKGFNKVVLKHLKKTKNLSSNEAISLEEIRIRVVQMLGSLGGQINKNLLTVTSSDEMMKSYVAWDREKRLSFAVPFREMKPVIFLDVFLPRVTELALTASDRQTKVAACELLHSMVMFMLGKATQMPEGGQGAPPMYQLYKRTFPVLLRLACDVDQVTRQLYEPLVMQLIHWFTNNKKFESQDTVALLEAILDGIVDPVDSTLRDFCGRCIREFLKWSIKQITPQQQEKSPVNTKSLFKRLYSLALHPNAFKRLGASLAFNNIYREFREEESLVEQFVFEALVIYMESLALAHADEKSLGTIQQCCDAIDHLCRIIEKKHVSLNKAKKRRLPRGFPPSASLCLLDLVKWLLAHCGRPQTECRHKSIELFYKFVPLLPGNRSPNLWLKDVLKEEGVSFLINTFEGGGCGQPSGILAQPTLLYLRGPFSLQATLCWLDLLLAALECYNTFIGERTVGALQVLGTEAQSSLLKAVAFFLESIAMHDIIAAEKCFGTGAAGNRTSPQEGERYNYSKCTVVVRIMEFTTTLLNTSPEGWKLLKKDLCNTHLMRVLVQTLCEPASIGFNIGDVQVMAHLPDVCVNLMKALKMSPYKDILETHLREKITAQSIEELCAVNLYGPDAQVDRSRLAAVVSACKQLHRAGLLHNILPSQSTDLHHSVGTELLSLVYKGIAPGDERQCLPSLDLSCKQLASGLLELAFAFGGLCERLVSLLLNPAVLSTASLGSSQGSVIHFSHGEYFYSLFSETINTELLKNLDLAVLELMQSSVDNTKMVSAVLNGMLDQSFRERANQKHQGLKLATTILQHWKKCDSWWAKDSPLETKMAVLALLAKILQIDSSVSFNTSHGSFPEVFTTYISLLADTKLDLHLKGQAVTLLPFFTSLTGGSLEELRRVLEQLIVAHFPMQSREFPPGTPRFNNYVDCMKKFLDALELSQSPMLLELMTEVLCREQQHVMEELFQSSFRRIARRGSCVTQVGLLESVYEMFRKDDPRLSFTRQSFVDRSLLTLLWHCSLDALREFFSTIVVDAIDVLKSRFTKLNESTFDTQITKKMGYYKILDVMYSRLPKDDVHAKESKINQVFHGSCITEGNELTKTLIKLCYDAFTENMAGENQLLERRRLYHCAAYNCAISVICCVFNELKFYQGFLFSEKPEKNLLIFENLIDLKRRYNFPVEVEVPMERKKKYIEIRKEAREAANGDSDGPSYMSSLSYLADSTLSEEMSQFDFSTGVQSYSYSSQDPRPATGRFRRREQRDPTVHDDVLELEMDELNRHECMAPLTALVKHMHRSLGPPQGEEDSVPRDLPSWMKFLHGKLGNPIVPLNIRLFLAKLVINTEEVFRPYAKHWLSPLLQLAASENNGGEGIHYMVVEIVATILSWTGLATPTGVPKDEVLANRLLNFLMKHVFHPKRAVFRHNLEIIKTLVECWKDCLSIPYRLIFEKFSGKDPNSKDNSVGIQLLGIVMANDLPPYDPQCGIQSSEYFQALVNNMSFVRYKEVYAAAAEVLGLILRYVMERKNILEESLCELVAKQLKQHQNTMEDKFIVCLNKVTKSFPPLADRFMNAVFFLLPKFHGVLKTLCLEVVLCRVEGMTELYFQLKSKDFVQVMRHRDDERQKVCLDIIYKMMPKLKPVELRELLNPVVEFVSHPSTTCREQMYNILMWIHDNYRDPESETDNDSQEIFKLAKDVLIQGLIDENPGLQLIIRNFWSHETRLPSNTLDRLLALNSLYSPKIEVHFLSLATNFLLEMTSMSPDYPNPMFEHPLSECEFQEYTIDSDWRFRSTVLTPMFVETQASQGTLQTRTQEGSLSARWPVAGQIRATQQQHDFTLTQTADGRSSFDWLTGSSTDPLVDHTSPSSDSLLFAHKRSERLQRAPLKSVGPDFGKKRLGLPGDEVDNKVKGAAGRTDLLRLRRRFMRDQEKLSLMYARKGVAEQKREKEIKSELKMKQDAQVVLYRSYRHGDLPDIQIKHSSLITPLQAVAQRDPIIAKQLFSSLFSGILKEMDKFKTLSEKNNITQKLLQDFNRFLNTTFSFFPPFVSCIQDISCQHAALLSLDPAAVSAGCLASLQQPVGIRLLEEALLRLLPAELPAKRVRGKARLPPDVLRWVELAKLYRSIGEYDVLRGIFTSEIGTKQITQSALLAEARSDYSEAAKQYDEALNKQDWVDGEPTEAEKDFWELASLDCYNHLAEWKSLEYCSTASIDSENPPDLNKIWSEPFYQETYLPYMIRSKLKLLLQGEADQSLLTFIDKAMHGELQKAILELHYSQELSLLYLLQDDVDRAKYYIQNGIQSFMQNYSSIDVLLHQSRLTKLQSVQALTEIQEFISFISKQGNLSSQVPLKRLLNTWTNRYPDAKMDPMNIWDDIITNRCFFLSKIEEKLTPLPEDNSMNVDQDGDPSDRMEVQEQEEDISSLIRSCKFSMKMKMIDSARKQNNFSLAMKLLKELHKESKTRDDWLVSWVQSYCRLSHCRSRSQGCSEQVLTVLKTVSLLDENNVSSYLSKNILAFRDQNILLGTTYRIIANALSSEPACLAEIEEDKARRILELSGSSSEDSEKVIAGLYQRAFQHLSEAVQAAEEEAQPPSWSCGPAAGVIDAYMTLADFCDQQLRKEEENASVIDSAELQAYPALVVEKMLKALKLNSNEARLKFPRLLQIIERYPEETLSLMTKEISSVPCWQFISWISHMVALLDKDQAVAVQHSVEEITDNYPQAIVYPFIISSESYSFKDTSTGHKNKEFVARIKSKLDQGGVIQDFINALDQLSNPELLFKDWSNDVRAELAKTPVNKKNIEKMYERMYAALGDPKAPGLGAFRRKFIQTFGKEFDKHFGKGGSKLLRMKLSDFNDITNMLLLKMNKDSKPPGNLKECSPWMSDFKVEFLRNELEIPGQYDGRGKPLPEYHVRIAGFDERVTVMASLRRPKRIIIRGHDEREHPFLVKGGEDLRQDQRVEQLFQVMNGILAQDSACSQRALQLRTYSVVPMTSRLGLIEWLENTVTLKDLLLNTMSQEEKAAYLSDPRAPPCEYKDWLTKMSGKHDVGAYMLMYKGANRTETVTSFRKRESKVPADLLKRAFVRMSTSPEAFLALRSHFASSHALICISHWILGIGDRHLNNFMVAMETGGVIGIDFGHAFGSATQFLPVPELMPFRLTRQFINLMLPMKETGLMYSIMVHALRAFRSDPGLLTNTMDVFVKEPSFDWKNFEQKMLKKGGSWIQEINVAEKNWYPRQKICYAKRKLAGANPAVITCDELLLGHEKAPAFRDYVAVARGSKDHNIRAQEPESGLSEETQVKCLMDQATDPNILGRTWEGWEPWM"

# Restraint Slopes
length_slope = 0.1
xl_slope = 0.1
semet_slope = 0.1
psipred_slope = 0.1


m = IMP.Model()
######################################
pdbfile = "../5luq_A_CA.pdb"
root_hier = IMP.atom.read_pdb(pdbfile, m)

# Define Structural Elements
# These are the residues in the PDB that correspond to structural elements.
# (start_res, length, SSID)
#-------------------------
elements=[(2602,14,'H'),(2623,25,'H'), (2655,10,'H')]

se = []

# Create list of StructureElements
for e in elements:
	se.append(setup_structural_element(root_hier, e))


#######################
# Set up Sequence
#######################
seq_chain = IMP.atom.Chain.setup_particle(IMP.Particle(m), "S")
root_hier.add_child(seq_chain)

res_particles = []
for i in range(len(seq)):
    pr = IMP.Particle(m)
    res = IMP.atom.Residue.setup_particle(pr,
                                        IMP.pmi.tools.get_residue_type_from_one_letter_code(seq[i]),
                                        i+1)
    res_particles.append(res.get_particle())
    IMP.atom.Mass.setup_particle(res.get_particle(), IMP.atom.get_mass(res.get_residue_type()))
    IMP.core.XYZR.setup_particle(res.get_particle())

    sel = IMP.atom.Selection(root_hier.get_children()[0], residue_index=i).get_selected_particles()
    if len(sel) == 1:
        IMP.core.XYZ(pr).set_coordinates(IMP.core.XYZ(sel[0]).get_coordinates())
        IMP.core.XYZ(pr).set_coordinates_are_optimized(True)

    seq_chain.add_child(res)

rests = []
#######################
# Set up Scoring Function
#######################

#xl_22 = [(2746,99,42),(2738,520,42),(2764,810,42),(2754,810,42),(2764,838,42),(2738,2003,42),(2694,2003,42),(2717,2107,42),(2738,2313,42),(2694,2445,42)]
xl_22 = [(2746,99,42),(2764,810,42),(2754,810,42),(2764,838,42),(2738,2003,42),(2694,2003,42),(2717,2107,42),(2738,2313,42),(2694,2445,42)]
xl_12 = [(2746,99,32),(2746,357,32),(2764,810,32),(2738,2227,32),(2746,2313,32)]
#xl_8 = [(2746,357,27),(2738,518,27),(2738,2001,27),(2717,2107,27),(2738,2227,27),(2746,2313,27)]
xl_8 = [(2746,357,27),(2738,2001,27),(2717,2107,27),(2738,2227,27),(2746,2313,27)]


import subprocess
xl_sites = []
for xl in xl_22 + xl_12 + xl_8:
#    xl_sites.append(xl[0])

    p = IMP.atom.Selection(root_hier, residue_index=xl[1]).get_selected_particles()
    if len(p) > 0:
        print IMP.core.XYZ(p[0])
        xl_sites.append(xl)

# XL Restraint
# ---------
# (res1, res2, length)
# Restraint is calculated using distance-per-residue library

xl_constant = 1.0 # # of SDs higher than mean
xl_slope = 0.05   # K-value of the harmonic
xlr = []
for xl in xl_sites:

    p1 = IMP.atom.Selection(root_hier, chain_id='S', residue_index=xl[0]).get_selected_particle_indexes()[0]
    p2 = IMP.atom.Selection(root_hier, chain_id='S', residue_index=xl[1]).get_selected_particle_indexes()[0]

    r = setup_conditional_pair_restraint(p1, p2, xl[2], xl_slope, xl_constant)
    pis = r.get_sequence_residue_particles()

    r1 = IMP.atom.Residue(m, pis[0])
    r2 = IMP.atom.Residue(m, pis[1])
    print "XL Setup between residues: ", xl[0], xl[1], " at a length of ", xl[2], "|", r1, r2

    rests.append(r)
    xlr.append(r)


# Completeness Restraint
# ---------
for s in se:
    r = setup_length_restraint(s, length_slope=length_slope)
    rests.append(r)


# Set up Terminal Elements
# For the loop-length distance restraint, we need to define SEs for the structured residues at
# each end of the unstructured domain(s)
terms = [[(27.365, -37.659, 9.827),2575],[(49.775, -42.147, 7.299),2774]]
term_ses = []
for t in terms:
    term_ses.append(setup_terminal_element(root_hier, t[0], t[1]))


# Loop Length distances
# ---------

# A structural domain sets the sequence of SEs.  Use the sequence to 
# add the E2E restraint.  model distance = last coord of SE1 - first coord of SE2.
# evaluated distance = (first resid of S2 - last resid of SE1) * res_dist

# This restraint is incompatible with SE swapping moves, as that changes the order of the SEs
# and, thus, would require re-formulating these restraints

# First, get the permutation from the 

print "THIS PERM: ", this_perm, perm_number
se_pairs = []
se_pairs.append((term_ses[0].get_particle_index(), se[this_perm[0]].get_particle_index()))
se_pairs.append((se[this_perm[0]].get_particle_index(), se[this_perm[1]].get_particle_index()))
se_pairs.append((se[this_perm[1]].get_particle_index(), se[this_perm[2]].get_particle_index()))
se_pairs.append((se[this_perm[2]].get_particle_index(), term_ses[-1].get_particle_index()))


secrs = []
for s in se_pairs:#[::-1]:
    #print s
    #print "> ", IMP.threading.StructureElement(m, s[0]).get_resindex_list()[-1], IMP.threading.StructureElement(m, s[0]).get_resindex_list()
    #print "> ", IMP.threading.StructureElement(m, s[1]).get_resindex_list()[0], IMP.threading.StructureElement(m, s[1]).get_resindex_list()
    
    r = add_SECR(s[0], s[1])
    secrs.append(r)
    #print ">>>> NRES:", r.get_number_of_residues(), r.get_model_distance(), r.get_max_distance()#, r.unprotected_evaluate(None), " "

    rests.append(r)


# PSIPRED Restraint
# ---------

# Read PSIPRED output
read_psipred_ss2("../data/psipred_dnapkcs.txt", res_particles)

# Get SE particle indexes
se_part_indexes = [s.get_particle_index() for s in se]

r = IMP.threading.SecondaryStructureParsimonyRestraint(m, se_part_indexes, seq_chain.get_particle_index(), 1.0)

rests.append(r)

sf = IMP.core.RestraintsScoringFunction(rests)


'''

# Add structural hierarchy
struct_hier = IMP.atom.Chain.setup_particle(IMP.Particle(m), "X")
IMP.core.XYZR.setup_particle(struct_hier)
IMP.atom.Mass.setup_particle(struct_hier, 1)
root_hier.add_child(struct_hier)

p = IMP.Particle(m)
f = IMP.atom.Fragment.setup_particle(p, range(6,13))

# The "correct" alignment results from Key values:
#  [(3,8,1,0), (19,19,1,0)]
'''

#########################
# Set up MC and Structure Elements Movers
#########################
mc = IMP.core.MonteCarlo(m)
mc.set_scoring_function(sf)
mc.set_return_best(False)

sems = []

for s in se:
    coordinates = s.get_coordinates()
    sem = IMP.threading.StructureElementMover(m, s.get_particle_index(), seq_chain.get_particle())
    sem.transform_coordinates()
    sems.append(sem)
    mc.add_mover(sem)

# Re-order SE lengths to correspond to this permutation
this_se_lengths = [se_lengths[this_perm[0]], se_lengths[this_perm[1]], se_lengths[this_perm[2]]]

srs = enumerate_start_resis_3( (2577,2774-this_se_lengths[-1]), this_se_lengths)

a0 = time.time()
print "System setup:", a0-init

print "INIT EVALUATE", sf.evaluate(False), [(r.get_name(), r.evaluate(False)) for r in rests]

all_output=[]

se_lengths=(14,25,10)

new_srs=[]

a = time.time()
print "Num models", len(srs), "enum setup time:", a-a0


first = n_thread * n_per_thread
last = first + n_per_thread

srs_for_us = srs[first:last]

for s in srs_for_us:
   new_srs.append((s[this_perm.index(0)], s[this_perm.index(1)],s[this_perm.index(2)]))


print "Nums", first, last, len(srs_for_us), srs_for_us[0], srs_for_us[1]
print "New Nums", first, last, len(srs_for_us), new_srs[0], new_srs[1]

f = open("enumerate_multi_"+str(n_thread)+"_"+str(perm_number)+"_"+str(n_per_thread)+".dat", "w")
out_dict = enumerate_start_resis_3_eval(new_srs, se, sems, sf, f)
a1 = time.time()

print "Final "+str(len(out_dict['models'])) +" models added"
print "Total time ", a1-a, "TPM=", (a1-a)/len(new_srs)


json.dump(out_dict, f) 
out_dict={}

IMP.atom.destroy(root_hier)
print "destroyed"








