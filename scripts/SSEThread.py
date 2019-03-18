import IMP
import IMP.atom
import IMP.core
import IMP.threading


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

def setup_conditional_pair_restraint(m, p1, p2, length, xl_slope, constant):
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
    m = s.get_model()
    uf = IMP.core.Linear(s.get_number_of_coordinates(), -1*length_slope)
    sf = IMP.core.AttributeSingletonScore(uf, IMP.FloatKey("length"))
    r = IMP.core.SingletonRestraint(m, sf, s.get_particle())
    return r

def add_SECR(m, p1, p2, slope=1, n_sds=1.0):
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
    m = se_list[0].get_model()
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

    rests = sf.create_restraints()

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
        new_mod['model'] = [s.get_all_key_values() for s in ses]
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