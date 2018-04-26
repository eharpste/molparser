import json
from concept_formation.continuous_value import ContinuousValue
from pprint import pprint

radical_dict = {
    0:"no_radical",
    1:"singlet",
    2:"doublet",
    3:"triplet"
}

charge_dict = {
    0:"outside_limits",
    1:"+3",
    2:"+2",
    3:"+1",
    4:"doublet_radical",
    5:"-1",
    6:"-2",
    7:"-3"
}

stereo_parity_dict = {
    0:"not_stereo",
    1:"odd",
    2:"even",
    3:"unmarked"
}

h_count_dict = {
    0:"H0",
    1:"H0",
    2:"H1",
    3:"H2",
    4:"H3",
    5:"H4"
}

bond_type_dict = {
    1:"Single",
    2:"Double",
    3:"Triple",
    4:"Aromatic",
    5:"Single_or_Double",
    6:"Single_or_Aromatic",
    7:"Double_or_Aromatic",
    8:"Any"
}

single_bond_stereo_dict = {
    0:"Not_stereo",
    1:"Up",
    4:"Either",
    6:"Down"
}

double_bond_stereo_dict = {
    0:"Use_coordinates",
    3:"Cis_or_trans"
}

def parse_counts_line(line):
    """
    Parses the counts line of a molecule and returns it asd a dictionary

    aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv
    aaa = number of atoms (current max 255)* 
    bbb = number of bonds (current max 255)* 
    lll = number of atom lists (max 30)* 
    fff = (obsolete)
    ccc = chiral flag: 0=not chiral, 1=chiral
    sss = number of stext entries 
    xxx = (obsolete)
    rrr = (obsolete)
    ppp = (obsolete)
    iii = (obsolete)
    mmm = number of lines of additional properties,
    vvvvv = version for the format
    """
    ret = {}
    ret["aaa"] = int(float(line[0:3]))
    ret["bbb"] = int(float(line[3:6]))
    ret["lll"] = int(float(line[6:9]))
    ret["ccc"] = int(float(line[12:15]))
    ret["sss"] = int(float(line[15:18]))
    ret["mmm"] = int(float(line[18:21]))
    ret["vvvvv"] = line[-5:]
    return ret

def parse_atom_line(line):
    """
    Parses a line from the atom block and returns it as a dictionary

    xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee
    [0:10] xxxxx.xxxx = x-coordinate
    [10:20] yyyyy.yyyy = y-coordinate
    [20:30] zzzzz.zzzz = z-coordinate
    [31:34] aaa = atomic symbol
    [34:36] dd = mass difference, i.e. difference from standard mass
    [36:39] ccc = charge 0 = uncharged or value other than these, 1 = +3, 2 = +2, 3 = +1, 4 = doublet radical, 5 = -1, 6 = -2, 7 = -3
    [39:42] sss = atom stereo parity 0 = not stereo, 1 = odd, 2 = even, 3 = either or unmarked stereo center
    [42:45] hhh = INGORED hydrogen count +1
    [45:48] bbb = IGNORED stereo care box
    [48:51] vvv = valence
    [51:54] HHH = IGNORED H0 designator
    [54:57] rrr = Not used
    [57:60] iii = Not used
    [60:63] mmm = IGNORED atom-atom mapping number 1 - number of atoms
    [63:66] nnn = IGNORED inversion/retention flag g 0 = property not applied 1 = configuration is inverted,2 = configuration is retained
    [66:69] eee = IGNORED 0 = property not applied, 1 = change on atom must be exactly as shown
    """
    ret = {}
    ret["xxx"] = float(line[0:10])
    ret["yyy"] = float(line[10:20])
    ret["zzz"] = float(line[20:30])
    ret["aaa"] = line[31:34].strip()
    ret["dd"] =  int(float(line[34:36]))
    ret["ccc"] = int(float(line[36:39]))
    ret["sss"] = int(float(line[39:42]))
    ret["hhh"] = int(float(line[42:45]))
    ret["bbb"] = int(float(line[45:48]))
    ret["vvv"] = int(float(line[48:51]))
    ret["HHH"] = int(float(line[51:54]))
    ret["mmm"] = int(float(line[60:63]))
    ret["nnn"] = int(float(line[63:66]))
    ret["eee"] = int(float(line[66:69]))
    return ret

def parse_bond_line(line):
    """
    Parses a line from a bondblock and turns it into a dict

    111222tttsssxxxrrrccc
    111 = number of atom 1
    222 = number of atom 2
    ttt = bond type
    sss = bond stereo
    xxx = not used
    rrr = bond topology
    ccc = reacting center status
    """
    ret = {}
    ret["111"] = int(float(line[0:3]))
    ret["222"] = int(float(line[3:6]))
    ret["ttt"] = int(float(line[6:9]))
    ret["sss"] = int(float(line[9:12]))
    ret["xxx"] = int(float(line[12:15]))
    ret["rrr"] = int(float(line[15:18]))
    ret["ccc"] = int(float(line[18:21]))
    return ret

def parse_data_name(line):
    """
    Parses the name of a data item line, which will be used as an attribute name
    """
    first = line.index("<")+1
    last = line.rindex(">")
    return line[first:last]

def apply_m_chg(line,mol):
    """
    Parses a CHG line from the property block.

    This will 0-out the charge on any atom that is not listed.
    0123456789
    M CHGnn8 aaa vvv ...
    aaa An atom number to alter the charge for
    vvv The ammount of charge on the target atom
    """
    if len(line) == 0:
        return mol
    for k in mol:
        if isinstance(k,str) and k.startswith('?atom'):
            mol[k]["charge"] = "0"
    line = line.rsplit(r'\s+')[3:]
    for i in range(0,len(line),2):
        aaa = int(line[i])
        vvv = line[i+1]
        if vvv[0] == "-":
            mol["?atom%04d"%aaa]["charge"] = vvv    
        else:
            mol["?atom%04d"%aaa]["charge"] = "+"+vvv
    return mol

def apply_m_rad(line,mol):
    """
    Parses a RAD line from the property block

    M RADnn8 aaa vvv ...
    aaa An atom number to alter the radical value of
    vvv A key into the radical_dict to set the radical value to
    """
    if len(line) == 0:
        return mol
    for k in mol:
        if isinstance(k,str) and k.startswith('?atom'):
            mol[k]["charge"] = "0"
            mol[k]["radical"] = "no radical"

    line = line.rsplit(r'\s+')[3:]
    for i in range(0,len(line),2):
        aaa = int(line[i])
        vvv = int(line[i+1])
        mol["?atom%04d"%aaa]["radical"] = radical_dict[vvv]
    return mol

def apply_m_iso(line,mol):
    """
    Parses and applies an ISO line from the property block

    M ISOnn8 aaa vvv ...
    aaa An atom numer to alter the mass_dif of
    vvv The mass_diff value to set
    """
    if len(line) == 0:
        return mol
    for k in mol:
        if isinstance(k,str) and k.startswith('?atom'):
            mol[k]["mass_diff"] = "0"

    line = line.rsplit(r'\s+')[3:]
    for i in range(0,len(line),2):
        aaa = int(line[i])
        vvv = int(line[i]+1)
        mol["?atom%04d"%aaa]["mass_diff"] = "+"+radical_dict[vvv]
    return mol

def merge_m_lines(lines):
    if len(lines) == 0:
        return ''
    line0 = lines[0]
    for l in lines[1:]:
        line0 += l[8:]
    return line0

def parse_mol(lines, verbose = False, data_items=False, tuple_realtions=False):
    """
    Parse the provided molfile and return a structured object representation
    that can be read by TRESTLE.
    """
    mol = {}
    num_atoms = 0
    num_bonds = 0
    num_lists = 0

    mol["name"] = lines[0]
    mol["software"] = lines[1].strip()
    if lines[2]:
        mol["comment"] = lines[2].strip()

    ### COUNTS LINE ###
    c_line = parse_counts_line(lines[3])
    num_atoms = c_line["aaa"]
    num_bonds = c_line["bbb"]
    num_lists = c_line["lll"]

    atom_dex = 4 + num_atoms
    bond_dex = atom_dex + num_bonds 

    if verbose:
        if c_line["ccc"] == 1:
            mol["chiral"] = True
        elif c_line["ccc"] == 0:
            mol["chiral"] = False

    mol["version"] = c_line["vvvvv"]

    ### ATOM BLOCK ###
    for l,line in enumerate(lines[4:atom_dex]):
        atom = {}
        a_line = parse_atom_line(line)
        atom["x"] = a_line["xxx"]
        atom["y"] = a_line["yyy"]
        atom["z"] = a_line["zzz"]

        atom["symbol"] = a_line["aaa"]
        #stringify mass difference or leave numeric?
        atom["mass_diff"] = str(a_line["dd"])
        atom["charge"] = charge_dict[a_line["ccc"]]
        atom["stereo_parity"] = stereo_parity_dict[a_line["sss"]]

        if a_line["vvv"] == 0:
            atom["valence"] = "no marking"
        if a_line["vvv"] > 0 and a_line["vvv"] < 15:
            atom["valence"] = str(a_line["vvv"])
        if a_line["vvv"] == 15:
            atom["valence"] = "zero valence"

        if verbose:
            atom["hydrogen_count"] = h_count_dict[a_line["hhh"]]


        mol["?atom%04d"%(l+1)] = atom

    ### BOND BLOCK ###
    for l,line in enumerate(lines[atom_dex:bond_dex]):
        bond = ['bond']
        b_line = parse_bond_line(line)

        bond.append(bond_type_dict[b_line["ttt"]])
        if bond_type_dict[b_line["ttt"]] == "Single":
            bond.append(single_bond_stereo_dict[b_line["sss"]])
        elif bond_type_dict[b_line["ttt"]] == "Double":
            bond.append(double_bond_stereo_dict[b_line["sss"]])
        bond.append("?atom%04d"%b_line["111"])
        bond.append("?atom%04d"%b_line["222"])

        if tuple_realtions:
            mol[tuple(bond)] = True
        else:
            mol['('+' '.join(bond)+')'] = True

    ### PROPERTIES BLOCK ###

    Ms = []
    for line in lines[bond_dex:]:
        if line == 'M  END':
            break
        if line[0] == "M":
            Ms.append(line)
        else:
            break
    
    if len(Ms) > 0:
        apply_m_chg(merge_m_lines([l for l in Ms if l.startswith('M  CHG')]),mol)
        apply_m_rad(merge_m_lines([l for l in Ms if l.startswith('M  RAD')]),mol)
        apply_m_iso(merge_m_lines([l for l in Ms if l.startswith('M  ISO')]),mol)

    ### DATA ITEMS ###
    if data_items:
        data_dex = bond_dex + len(Ms)
        opened = False
        data_header = ""
        data_list = []
        for line in lines[data_dex:]:
            if line:
                if not opened and line[0] == ">":
                    opened = True
                    data_header = parse_data_name(line)
                    data_list = []
                    continue
                elif opened:
                    data_list.append(line)
            elif opened:
                mol[data_header] = '\n'.join(data_list)
                opened = False

    return mol

def parse_sdf_file(filename,verbose = False,data_items=False,tuple_realtions=False,n=-1):
    ret = []
    curr = []
    count = 0
    with open(filename,"r") as molefile:
        for line in molefile:
            line = line.rstrip('\r\n')
            if not line == "$$$$":
                curr.append(line)
            else:
                if len(curr) > 0:
                    ret.append(parse_mol(curr,verbose,data_items,tuple_realtions))    
                    count += 1
                    if n > 0 and count > n:
                        break
                curr = []
    return ret

if __name__ == '__main__':
    # mols = parse_sdf_file("./example.mol",True,True,True,n=100)
    mols = parse_sdf_file("./testSDFs/Compound_000000001_000025000.sdf",True,True,True,n=25000)
    pprint(mols[0])