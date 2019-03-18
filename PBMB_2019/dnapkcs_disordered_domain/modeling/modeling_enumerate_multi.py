from __future__ import print_function

import IMP
import numpy
import IMP.pmi.tools
import IMP.atom
import IMP.core
import IMP.threading
import json

import sys
sys.path.insert(0, '../../../scripts/')
from SSEThread import *

# These arguments facilitate the multi-threading of enumeration runs
# They divide the total list of models into n_per_thread chunks and assigns
# chunk n_thread to this instance.
n_per_thread = int(sys.argv[1])
n_thread = int(sys.argv[2])
perm_number = int(sys.argv[3]) # Permutation of the SE sequence

# Hard-coded in permutations for a sequence of three
perms = [(0,1,2), (0,2,1), (1,0,2), (1,2,0), (2,1,0), (2,0,1)]
this_perm = perms[perm_number]

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

    r = setup_conditional_pair_restraint(m, p1, p2, xl[2], xl_slope, xl_constant)
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
    r = add_SECR(m, s[0], s[1])
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

print("Total Score on initialization:", sf.evaluate(False), [(r.get_name(), r.evaluate(False)) for r in rests])

all_output=[]


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


