import fileinput
from shutil import copyfile
from collections import defaultdict

'''
This script takes an SSEThread output and produces a PDB for each model in that output

Using the pdb_template file, it grabs lines from the PDB corresponding to the 
structure elements (elements).  Using the SSEThread output (input_file),
a new line with the correct residue assignment is created.

A new PDB is created swapping this new line for the old line for all residues in the
structure elements. Requires that the pdb_template has 'UNK' for all element residue names.

The new PDB is written to the file indicated by output_prefix, plus the model number. 
'''


input_file = "bsms_10000.dat"
pdb_template = "../data/5luq_A_CA.pdb"
elements = [(2602,14,'H'),(2623,25,'H'), (2655,10,'H')]
outpdb_prefix = "./bsms_pdbs/bsms_"


# The following class is copied from IMP.pmi.tools on 8/13/2019
# https://github.com/salilab/imp/blob/758929befa2dfa453008a958b439972805ddaec3/modules/pmi/pyext/src/tools.py

class ThreeToOneConverter(defaultdict):
    """This class converts three to one letter codes, and return X for any unknown codes"""
    def __init__(self,is_nucleic=False):

        if not is_nucleic:
            threetoone = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
                              'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
                              'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
                              'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
                              'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V', 'UNK': 'X'}
        else:
            threetoone = {'ADE': 'A', 'URA': 'U', 'CYT': 'C', 'GUA': 'G',
                              'THY': 'T', 'UNK': 'X'}

        defaultdict.__init__(self,lambda: "X", threetoone)

# The following function is copied from IMP.pmi.tools on 8/13/2019
# https://github.com/salilab/imp/blob/758929befa2dfa453008a958b439972805ddaec3/modules/pmi/pyext/src/tools.py

def get_residue_type_from_one_letter_code(code,is_nucleic=False):
    threetoone=ThreeToOneConverter(is_nucleic)
    one_to_three={}
    for k in threetoone:
        one_to_three[threetoone[k]] = k
    return one_to_three[code]


def replace_unk_with_res_code_and_index(text_line, residue, old_resnum, new_resnum):
    '''
    Replaces all instances of "UNK" in the given string 
    with the 3 letter residue code of the given amino acid letter
    '''
    res_code = get_residue_type_from_one_letter_code(residue)
    new_text_line = text_line.replace("UNK", res_code)
    new_text_line = new_text_line.replace(str(old_resnum), str(new_resnum))

    return new_text_line


def get_residue_line_from_pdb(file_handle, residue_number, chainID="A", atom_name="CA"):

    with open(file_handle) as search:
        for lnie in search:
            lnie = lnie.rstrip()  # remove '\n' at end of line

            #print(lnie[22:26].strip(), "||", residue_number, atom_name, chainID, "||", lnie[12:16].strip(),lnie[21], lnie)

            if lnie[0:6]=="ATOM  ":
                if int(lnie[22:26].strip()) == residue_number:
                    #print("JIDOKSNO", lnie[12:16].strip(), atom_name, lnie[12:16].strip() == atom_name)
                    #exit()
                    if lnie[12:16].strip() == atom_name and lnie[21] == chainID:  
                        search.close()
                        return lnie
        search.close()
        print("No match found")
        return None



pdb_element_start_resis = ()

sequence = "MAGSGAGVRCSLLRLQETLSAADRCGAALAGHQLIRGLGQECVLSSSPAVLALQTSLVFSRDFGLLVFVRKSLNSIEFRECREEILKFLCIFLEKMGQKIAPYSVEIKNTCTSVYTKDRAAKCKIPALDLLIKLLQTFRSSRLMDEFKIGELFSKFYGELALKKKIPDTVLEKVYELLGLLGEVHPSEMINNAENLFRAFLGELKTQMTSAVREPKLPVLAGCLKGLSSLLCNFTKSMEEDPQTSREIFNFVLKAIRPQIDLKRYAVPSAGLRLFALHASQFSTCLLDNYVSLFEVLLKWCAHTNVELKKAALSALESFLKQVSNMVAKNAEMHKNKLQYFMEQFYGIIRNVDSNNKELSIAIRGYGLFAGPCKVINAKDVDFMYVELIQRCKQMFLTQTDTGDDRVYQMPSFLQSVASVLLYLDTVPEVYTPVLEHLVVMQIDSFPQYSPKMQLVCCRAIVKVFLALAAKGPVLRNCISTVVHQGLIRICSKPVVLPKGPESESEDHRASGEVRTGKWKVPTYKDYVDLFRHLLSSDQMMDSILADEAFFSVNSSSESLNHLLYDEFVKSVLKIVEKLDLTLEIQTVGEQENGDEAPGVWMIPTSDPAANLHPAKPKDFSAFINLVEFCREILPEKQAEFFEPWVYSFSYELILQSTRLPLISGFYKLLSITVRNAKKIKYFEGVSPKSLKHSPEDPEKYSCFALFVKFGKEVAVKMKQYKDELLASCLTFLLSLPHNIIELDVRAYVPALQMAFKLGLSYTPLAEVGLNALEEWSIYIDRHVMQPYYKDILPCLDGYLKTSALSDETKNNWEVSALSRAAQKGFNKVVLKHLKKTKNLSSNEAISLEEIRIRVVQMLGSLGGQINKNLLTVTSSDEMMKSYVAWDREKRLSFAVPFREMKPVIFLDVFLPRVTELALTASDRQTKVAACELLHSMVMFMLGKATQMPEGGQGAPPMYQLYKRTFPVLLRLACDVDQVTRQLYEPLVMQLIHWFTNNKKFESQDTVALLEAILDGIVDPVDSTLRDFCGRCIREFLKWSIKQITPQQQEKSPVNTKSLFKRLYSLALHPNAFKRLGASLAFNNIYREFREEESLVEQFVFEALVIYMESLALAHADEKSLGTIQQCCDAIDHLCRIIEKKHVSLNKAKKRRLPRGFPPSASLCLLDLVKWLLAHCGRPQTECRHKSIELFYKFVPLLPGNRSPNLWLKDVLKEEGVSFLINTFEGGGCGQPSGILAQPTLLYLRGPFSLQATLCWLDLLLAALECYNTFIGERTVGALQVLGTEAQSSLLKAVAFFLESIAMHDIIAAEKCFGTGAAGNRTSPQEGERYNYSKCTVVVRIMEFTTTLLNTSPEGWKLLKKDLCNTHLMRVLVQTLCEPASIGFNIGDVQVMAHLPDVCVNLMKALKMSPYKDILETHLREKITAQSIEELCAVNLYGPDAQVDRSRLAAVVSACKQLHRAGLLHNILPSQSTDLHHSVGTELLSLVYKGIAPGDERQCLPSLDLSCKQLASGLLELAFAFGGLCERLVSLLLNPAVLSTASLGSSQGSVIHFSHGEYFYSLFSETINTELLKNLDLAVLELMQSSVDNTKMVSAVLNGMLDQSFRERANQKHQGLKLATTILQHWKKCDSWWAKDSPLETKMAVLALLAKILQIDSSVSFNTSHGSFPEVFTTYISLLADTKLDLHLKGQAVTLLPFFTSLTGGSLEELRRVLEQLIVAHFPMQSREFPPGTPRFNNYVDCMKKFLDALELSQSPMLLELMTEVLCREQQHVMEELFQSSFRRIARRGSCVTQVGLLESVYEMFRKDDPRLSFTRQSFVDRSLLTLLWHCSLDALREFFSTIVVDAIDVLKSRFTKLNESTFDTQITKKMGYYKILDVMYSRLPKDDVHAKESKINQVFHGSCITEGNELTKTLIKLCYDAFTENMAGENQLLERRRLYHCAAYNCAISVICCVFNELKFYQGFLFSEKPEKNLLIFENLIDLKRRYNFPVEVEVPMERKKKYIEIRKEAREAANGDSDGPSYMSSLSYLADSTLSEEMSQFDFSTGVQSYSYSSQDPRPATGRFRRREQRDPTVHDDVLELEMDELNRHECMAPLTALVKHMHRSLGPPQGEEDSVPRDLPSWMKFLHGKLGNPIVPLNIRLFLAKLVINTEEVFRPYAKHWLSPLLQLAASENNGGEGIHYMVVEIVATILSWTGLATPTGVPKDEVLANRLLNFLMKHVFHPKRAVFRHNLEIIKTLVECWKDCLSIPYRLIFEKFSGKDPNSKDNSVGIQLLGIVMANDLPPYDPQCGIQSSEYFQALVNNMSFVRYKEVYAAAAEVLGLILRYVMERKNILEESLCELVAKQLKQHQNTMEDKFIVCLNKVTKSFPPLADRFMNAVFFLLPKFHGVLKTLCLEVVLCRVEGMTELYFQLKSKDFVQVMRHRDDERQKVCLDIIYKMMPKLKPVELRELLNPVVEFVSHPSTTCREQMYNILMWIHDNYRDPESETDNDSQEIFKLAKDVLIQGLIDENPGLQLIIRNFWSHETRLPSNTLDRLLALNSLYSPKIEVHFLSLATNFLLEMTSMSPDYPNPMFEHPLSECEFQEYTIDSDWRFRSTVLTPMFVETQASQGTLQTRTQEGSLSARWPVAGQIRATQQQHDFTLTQTADGRSSFDWLTGSSTDPLVDHTSPSSDSLLFAHKRSERLQRAPLKSVGPDFGKKRLGLPGDEVDNKVKGAAGRTDLLRLRRRFMRDQEKLSLMYARKGVAEQKREKEIKSELKMKQDAQVVLYRSYRHGDLPDIQIKHSSLITPLQAVAQRDPIIAKQLFSSLFSGILKEMDKFKTLSEKNNITQKLLQDFNRFLNTTFSFFPPFVSCIQDISCQHAALLSLDPAAVSAGCLASLQQPVGIRLLEEALLRLLPAELPAKRVRGKARLPPDVLRWVELAKLYRSIGEYDVLRGIFTSEIGTKQITQSALLAEARSDYSEAAKQYDEALNKQDWVDGEPTEAEKDFWELASLDCYNHLAEWKSLEYCSTASIDSENPPDLNKIWSEPFYQETYLPYMIRSKLKLLLQGEADQSLLTFIDKAMHGELQKAILELHYSQELSLLYLLQDDVDRAKYYIQNGIQSFMQNYSSIDVLLHQSRLTKLQSVQALTEIQEFISFISKQGNLSSQVPLKRLLNTWTNRYPDAKMDPMNIWDDIITNRCFFLSKIEEKLTPLPEDNSMNVDQDGDPSDRMEVQEQEEDISSLIRSCKFSMKMKMIDSARKQNNFSLAMKLLKELHKESKTRDDWLVSWVQSYCRLSHCRSRSQGCSEQVLTVLKTVSLLDENNVSSYLSKNILAFRDQNILLGTTYRIIANALSSEPACLAEIEEDKARRILELSGSSSEDSEKVIAGLYQRAFQHLSEAVQAAEEEAQPPSWSCGPAAGVIDAYMTLADFCDQQLRKEEENASVIDSAELQAYPALVVEKMLKALKLNSNEARLKFPRLLQIIERYPEETLSLMTKEISSVPCWQFISWISHMVALLDKDQAVAVQHSVEEITDNYPQAIVYPFIISSESYSFKDTSTGHKNKEFVARIKSKLDQGGVIQDFINALDQLSNPELLFKDWSNDVRAELAKTPVNKKNIEKMYERMYAALGDPKAPGLGAFRRKFIQTFGKEFDKHFGKGGSKLLRMKLSDFNDITNMLLLKMNKDSKPPGNLKECSPWMSDFKVEFLRNELEIPGQYDGRGKPLPEYHVRIAGFDERVTVMASLRRPKRIIIRGHDEREHPFLVKGGEDLRQDQRVEQLFQVMNGILAQDSACSQRALQLRTYSVVPMTSRLGLIEWLENTVTLKDLLLNTMSQEEKAAYLSDPRAPPCEYKDWLTKMSGKHDVGAYMLMYKGANRTETVTSFRKRESKVPADLLKRAFVRMSTSPEAFLALRSHFASSHALICISHWILGIGDRHLNNFMVAMETGGVIGIDFGHAFGSATQFLPVPELMPFRLTRQFINLMLPMKETGLMYSIMVHALRAFRSDPGLLTNTMDVFVKEPSFDWKNFEQKMLKKGGSWIQEINVAEKNWYPRQKICYAKRKLAGANPAVITCDELLLGHEKAPAFRDYVAVARGSKDHNIRAQEPESGLSEETQVKCLMDQATDPNILGRTWEGWEPWM"

# 1) Open bsms file
f = open(input_file, "r")

k = 0
# For each line in bsms, extract the model
for l in f.readlines():
    model = l.split("|")[0].strip()
    score = model.split(" ")[3]
    mod_start_resis = model.split(" ")[0:3]

    if len(mod_start_resis) != len(elements):
        raise Exception("Number of elements," +str(len(mod_start_resis))+", in bsms.dat is not the same as defined in the elements list of the make_pdbs.py script")

    # Create the line swapping dictionary
    swap_dict = {}

    for i in range(len(mod_start_resis)):

        # Start residue of the model SE
        sr = int(float(mod_start_resis[i]))

        # The starting residue and length of the element in the reference PDB
        elem_sr = int(elements[i][0])
        elem_len = int(elements[i][1])

        for j in range(elem_len):
            # Get sequence residue number
            res_number = j + elem_sr
            new_res_num = sr + j

            new_res_letter = sequence[res_number-1]

            # Get the residue coordinates from the PDB
            res_line = get_residue_line_from_pdb(pdb_template, res_number).rstrip('\r\n')

            # replace the UNK with the correct residue. This is the line we will add to the BSM pdb
            new_res_line = replace_unk_with_res_code_and_index(res_line, new_res_letter, res_number, new_res_num)+"  "

            swap_dict[res_line+"  "] = new_res_line

    # Copy the template file and swap line in the new file
    new_file = outpdb_prefix+str(k)+".pdb"

    print(new_file)#, list(swap_dict.keys())[0][5:38], "||||", swap_dict[list(swap_dict.keys())[0]][5:38])

    of = open(new_file, "w")
    inpdb = open(pdb_template, "r")

    for lnie in inpdb.readlines():
        
        if lnie.rstrip('\n\r') in swap_dict:
            of.write(swap_dict[lnie.rstrip('\n\r')]+'\n')
        elif "UNK" not in lnie:
            of.write(lnie)

    of.close()
    inpdb.close()

    k+=1

