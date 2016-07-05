import json
from pprint import pprint

radical_dict = {
    0:"no radical",
    1:"singlet",
    2:"doublet",
    3:"triplet"
}

charge_dict = {
    0:"outside limits",
    1:"+3",
    2:"+2",
    3:"+1",
    4:"doublet radical",
    5:"-1",
    6:"-2",
    7:"-3"
}

stereo_parity_dict = {
    0:"not stereo",
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
    5:"Single or Double",
    6:"Single or Aromatic",
    7:"Double or Aromatic",
    8:"Any"
}

single_bond_stereo_dict = {
    0:"Not stereo",
    1:"Up",
    4:"Either",
    6:"Down"
}

double_bond_stereo_dict = {
    0:"Use coordinates",
    3:"Cis or trans"
}

def parse_counts_line(line):
    #aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv
    #aaa = number of atoms (current max 255)* 
    #bbb = number of bonds (current max 255)* 
    #lll = number of atom lists (max 30)* 
    #fff = (obsolete)
    #ccc = chiral flag: 0=not chiral, 1=chiral
    #sss = number of stext entries 
    #xxx = (obsolete)
    #rrr = (obsolete)
    #ppp = (obsolete)
    #iii = (obsolete)
    #mmm = number of lines of additional properties,
    #vvvvv = version for the format
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
    #xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee
    #[0:10] xxxxx.xxxx = x-coordinate
    #[10:20] yyyyy.yyyy = y-coordinate
    #[20:30] zzzzz.zzzz = z-coordinate
    #[31:34] aaa = atomic symbol
    #[34:36] dd = mass difference, i.e. difference from standard mass
    #[36:39] ccc = charge 0 = uncharged or value other than these, 1 = +3, 2 = +2, 3 = +1, 4 = doublet radical, 5 = -1, 6 = -2, 7 = -3
    #[39:42] sss = atom stereo parity 0 = not stereo, 1 = odd, 2 = even, 3 = either or unmarked stereo center
    #[42:45] hhh = INGORED hydrogen count +1
    #[45:48] bbb = IGNORED stereo care box
    #[48:51] vvv = valence
    #[51:54] HHH = IGNORED H0 designator
    #[54:57] rrr = Not used
    #[57:60] iii = Not used
    #[60:63] mmm = IGNORED atom-atom mapping number 1 - number of atoms
    #[63:66] nnn = IGNORED inversion/retention flag g 0 = property not applied 1 = configuration is inverted,2 = configuration is retained
    #[66:69] eee = IGNORED 0 = property not applied, 1 = change on atom must be exactly as shown
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
    #111222tttsssxxxrrrccc
    #111 = number of atom 1
    #222 = number of atom 2
    #ttt = bond type
    #sss = bond stereo
    #xxx = not used
    #rrr = bond topology
    #ccc = reacting center status
    ret = {}
    ret["111"] = int(float(line[0:3]))
    ret["222"] = int(float(line[3:6]))
    ret["ttt"] = int(float(line[6:9]))
    ret["sss"] = int(float(line[9:12]))
    ret["xxx"] = int(float(line[12:15]))
    ret["rrr"] = int(float(line[15:18]))
    ret["ccc"] = int(float(line[18:21]))
    return ret

def parse_prop_name(line):
    first = line.index("<")+1
    last = line.rindex(">")
    return line[first:last]

def apply_m_chg(line,mol):
    #print('apply_m_chg')
    for k in mol:
        if k.startswith('?atom'):
        #if k[0:5] == "?atom":
            mol[k]["charge"] = "0"
    #print(line)
    #print(line[9:])
    for i in range(9,len(line),8):
        aaa = line[i:i+4].strip()
        vvv = line[i+4:i+8].strip()
        #print(str(i),aaa,vvv)
        #pprint(mol)
        if vvv[0] == "-":
            mol["?atom"+aaa]["charge"] = vvv    
        else:
            mol["?atom"+aaa]["charge"] = "+"+vvv
#123456789A123456789B123456789
#M  CHG  2   2  -1   5   1
def apply_m_rad(line,mol):
    for k in mol:
        if k.startswith('?atom'):
        #if k[0:4] == "atom":
            mol[k]["charge"] = "0"
            mol[k]["radical"] = "no radical"

    for i in range(11,len(line),6):
        aaa = line[i:i+3].strip()
        vvv = int(float(line[i+3:i+6].strip()))
        mol["?atom"+aaa]["radical"] = radical_dict[vvv]

def apply_m_iso(line,mol):
    for k in mol:
        if k.startswith('?atom'):
        #if k[0:4] == "atom":
            mol[k]["mass_diff"] = "0"

    for i in range(11,len(line),6):
        aaa = line[i:i+3].strip()
        vvv = line[i+3:i+6].strip()
        mol["?atom"+aaa]["mass_diff"] = "+"+radical_dict[vvv]


def parse_mol(lines, verbose = False, properties=False):
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


        mol["?atom"+str(l+1)] = atom

    for l,line in enumerate(lines[atom_dex:bond_dex]):
        bond = ['bond']
        b_line = parse_bond_line(line)

        bond.append(bond_type_dict[b_line["ttt"]])
        if bond_type_dict[b_line["ttt"]] == "Single":
            bond.append(single_bond_stereo_dict[b_line["sss"]].replace(' ','_'))
        elif bond_type_dict[b_line["ttt"]] == "Double":
            bond.append(double_bond_stereo_dict[b_line["sss"]].replace(' ','_'))
        bond.append("?atom"+str(b_line["111"]))
        bond.append("?atom"+str(b_line["222"]))

        mol['('+' '.join(bond)+')'] = True

    Ms = []
    for line in lines[bond_dex:]:
        if line[0] == "M":
            Ms.append(line)
        else:
            break
    
    if len(Ms) > 1:
        if len(Ms) > 2:
            for m in Ms[1:-1]:
                Ms[0] += m[11:]
        
        if Ms[0][3:6] == "CHG":
            apply_m_chg(Ms[0],mol)

        elif Ms[0][3:6] == "RAD":
            apply_m_rad(Ms[0],mol)

        elif Ms[0][3:6] == "ISO":
            apply_m_iso(Ms[0],mol)

    if properties:
        prop_dex = bond_dex + len(Ms)
        opened = False
        prop_header = ""
        prop_list = []
        for line in lines[prop_dex:]:
            if line:
                if not opened and line[0] == ">":
                    opened = True
                    prop_header = parse_prop_name(line)
                    prop_list = []
                    continue
                elif opened:
                    prop_list.append(line)
            elif opened:
                mol[prop_header] = '\n'.join(prop_list)
                opened = False

    return mol

def parse_molefile(filename,verbose = False,properties=False,n=-1):
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
                    ret.append(parse_mol(curr,verbose,properties))    
                    count += 1
                    if n > 0 and count > n:
                        break
                curr = []
    return ret

if __name__ == '__main__':
    mols = parse_molefile("./example.mol",True,True,n=5)
    #mols = parse_molefile("./testSDFs/Compound_000000001_000025000.sdf",True,True,n=1)
    pprint(mols)