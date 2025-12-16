"""
# -----------------------------------------------------------------------------
# New/updated version of XDStrategy by Frej (original by Lennard Krause). 
# 
# Notes:
# Remember to update the location of XDConstraints if this is to be used
# directly through this script.
# 
# Latest changes:
# - Added version number :D
# - Moved path to XDConstraints for easy access
# 
# Todo (not a prioritised list):
# - Set *model correctly (-1 (or -2?) for IAM)
# - Incorporate refining κ on either core or valence
# - Make standard sequence actually usable
# - Use standard sequence directly
# - Make sure to always save the original master file
# - Make SITESYM actually do something:
#   - Geneate symmetry constraints on XYZ and U
#   - Generate pseudosymm for MPs
# - Make some functions and remove some indents:
#   - Scalefactor handling
# - Actually make the script decent (and borderline PEP8 compliant)
# - Write some documentation/docstrings (perhaps)
# -----------------------------------------------------------------------------
"""

import re 
import sys
import os
from subprocess import run
from copy import deepcopy
from glob import glob
from pathlib import Path
from collections import defaultdict

__version__ = 'v0.0.2, 16.12.2025'

standard = [
'# Instructions for masterfiles',
'# Masterfiles will be numbered automatically xd01.mas, ...',
'# Specifications are made with [] and seperated by ; (eg. U3[Co;Br])'
'# Exclude atoms with ! (eg. !H to exclude hydrogens)',
'# In- and exclusions can be element or atom specific (eg. Co or Co(1))'
'# Inclusion and exclusion are mutually exclusive (only use one!)',
'# XYZ and U2-4 are ALWAYS tagged with [!H] (use HXYZ or HU1 instead)'
'# Possible instructions (WIP):',
'# XYZ, HXYZ, U1, U2, U3, U4, HU1, M, D, Q, O, H, CC,',
'# SCALE, SINTHL, SIGOBS, KAPPA, KAPPAP, CON',
'SCALE',
'SCALE CC M',
'SCALE CC M D Q O[!H] H[!H]',
'SCALE CC M D Q O[!H] H[!H] U2',
'SCALE CC M D Q O[!H] H[!H] XYZ U2',
'SCALE CC M D Q O[!H] H[!H] XYZ U2 KAPPA',
'SCALE CC HXYZ',
'SCALE CC M D Q O[!H] H[!H] XYZ U2 KAPPA',
'SCALE CC KAPPAP',
'SCALE CC M D Q O[!H] H[!H] XYZ U2 KAPPA',
'SCALE CC SIGOBS[0] M D Q O[!H] H[!H] XYZ U2 KAPPA'
]

# Path to XDConstraints script (by Lennard Krause)
XDCon_path = r''

## Patterns
atom_table_regex = re.compile(r'(?P<ATOM>[a-zA-Z\(\)0-9]+)\s+(?P<ATOM0>[a-zA-Z\(\)0-9]+)\s+(?P<AX1>[XYZxyz])\s+(?P<ATOM1>[a-zA-Z\(\)0-9]+)\s+(?P<ATOM2>[a-zA-Z\(\)0-9]+)\s+(?P<AX2>[XYZxyz])\s+(?P<RL>[RLrl])\s+(?P<TP>\d)\s+(?P<TBL>\d+)\s+(?P<KAP>\d+)\s+(?P<LMX>\d)\s+(?P<SITESYM>[0-9a-zA-Z_]*){0,1}(\s+(?P<CHEMCON>.*)\s*)*\n')
key_table_regex  = re.compile(r'(?P<ATOM>[a-zA-Z\(\)0-9]+)\s+(?P<XYZ>[01]{3})\s+(?P<U2>[01]{6})\s+(?P<U3>[01]{10})\s+(?P<U4>[01]{15})\s+(?P<M>[01]{2})\s+(?P<D>[01]{3})\s+(?P<Q>[01]{5})\s+(?P<O>[01]{7})\s+(?P<H>[01]{9}).*\n')
scat_table_regex = re.compile(r'(?P<SCAT>[a-zA-Z\(\)0-9\+\-]+)\s+(?P<CORE>\w+)\s+(?P<SPHV>\w+)\s+(?P<DEFV>\w+)\s+(?P<S1>[0-9\-]+)\s+(?P<S2>[0-9\-]+)\s+(?P<S3>[0-9\-]+)\s+(?P<S4>[0-9\-]+)\s+(?P<P2>[0-9\-]+)\s+(?P<P3>[0-9\-]+)\s+(?P<P4>[0-9\-]+)\s+(?P<D3>[0-9\-]+)\s+(?P<D4>[0-9\-]+)\s+(?P<F4>[0-9\-]+)\s+(?P<S5>[0-9\-]+)\s+(?P<P5>[0-9\-]+)\s+(?P<S6>[0-9\-]+)\s+(?P<P6>[0-9\-]+)\s+(?P<D5>[0-9\-]+)\s+(?P<S7>[0-9\-]+)\s+(?P<D6>[0-9\-]+)\s+(?P<F5>[0-9\-]+)\s+(?P<DELF>[0-9\.\-]+)\s+(?P<DELFf>[0-9\.\-]+)\s+(?P<NSCTL>[0-9\.\-]+)\s*')
kappa_regex      = re.compile(r'KAPPA\s+[01]{6}\s*')
four_regex       = re.compile(r'!{0,1}FOUR\s+fmod1\s+-?\d\s+\d\s+\d\s+\d\s+fmod2\s+-?\d\s+\d\s+\d\s+\d')
model_regex      = re.compile(r'SELECT\s+\**model\s+4\s+[23]\s+1\s+0\s+based_on\s+F\^2\s+test\s+verbose\s+1\s+')
skip_regex       = re.compile(r'!{0,1}SKIP\s+[\*]{0,1}obs\s+[d\d.]+\s+[d\d.]+\s+[\*]{0,1}sigobs\s+[d\d.]+\s+[d\d.]+\s+[\*]{0,1}sinthl\s+[d\d.]+\s+[d\d.]+\s+')
vibcon_regex     = re.compile(r'CON(\s+[\d.-]+\s+[U\d/]+)+\s+=\s+0\s+.*')
reset_regex      = re.compile(r'RESET\s+BOND(\s+[a-zA-Z\(\)0-9]+){2}\s+[\d.]+\s+')
cycles_regex     = re.compile(r'SELECT\s+cycle\s+[0-9\-]+\s+dampk\s+[0-9\.]+\s+cmin\s+[0-9\.]+\s+cmax\s+[0-9\.]+\s+eigcut\s+[d0-9\.\-]+\s+\*?convcrit\s+[d0-9\.\-]+\s+')
scale_regex      = re.compile(r'^SCALE\s+\d+')
dum_regex        = re.compile(r'^DUM\d+\s+-?\d+.\d+\s+-?\d+.\d+\s+-?\d+.\d+')

## Find files
def get_files():
    print()  # Beauty spacing
    # Masterfile
    while True:
        master_input = input('Enter name of source XD master file [xd.mas]: ') or 'xd.mas'
        try:
            oMasName = master_input.split('.')[0] + '.mas'
            with open(oMasName, 'r') as ifile:
                oMas = ifile.readlines()
            break
        except FileNotFoundError: print(f'Source XD masterfile ({oMasName}) not found. Try again.')

    # Instructions file
    inst_file_input = input('Enter name of instructions file (.inst) [xd.inst]: ') or 'xd.inst'
    try:
        with open(inst_file_input.split('.')[0] + '.inst', 'r') as ifile:
            inst_file = ifile.readlines()
    except FileNotFoundError:
        WriteInst = input('Instructions file not found. Write standard? (Y/N) [Y]: ') or 'Y'
        if WriteInst.upper() == 'Y':
            with open('xd.inst', 'w') as TempFile:
                for i in standard: TempFile.write(i + '\n')
            print('Instructions written to xd.inst.')
            Cont = (input('Set up masterfiles with standard strategy? (Y/N) [N]: ') or 'N').upper()
            if Cont == 'Y':
                print('You forgot to add this feature.')
                inst_file = standard
            else:
                print('Please restart program with the correct instructions file.')
                sys.exit()

        else:
            print('Please restart program with an instructions file.')
            sys.exit()

    # Pruning comments and empty lines in instructions:
    inst_file = [line for line in inst_file if not (line.startswith('#') or re.match(r'\s*\n', line))]

    # Constraints
    Consts = input('Load (Y/N) or Write (W) constraints file(s)? [Y]: ') or 'Y'
    print()  # Beauty spacing
    
    if Consts.upper() == 'Y':
        const_files = []
        temp_const_files = glob(r'*.con') + glob(r'*.const')
        load_all = input(f"Found these constraint files: [{', '.join(temp_const_files)}]. Should they all be loaded? (Y/N) [Y]: ") or 'Y'
        if load_all.upper() == 'Y':
            const_files = temp_const_files

        else:
            while True:
                const_files_input = input('Input constraint file(s). Separate files with space: ')
                if const_files_input:
                    files = [file for file in const_files_input.split(' ') if file]

                    try:
                        for file in files:
                            with open(file, 'r'):
                                const_files.append(file)
                                print(f'{file} loaded.')
                        break
                    except FileNotFoundError:
                        print(f'{file} not found.')
                else:
                    print('Skipped loading constraints.')
                    break

    elif Consts.upper() == 'W' and Path(XDCon_path).exists():
        existing_cons = glob(r'*.con') + glob(r'*.const')
        print(f'Running XDConstraints on {oMasName} to write constraints file.')
        XDCon = run(
                ['python', XDCon_path],
                input=oMasName,
                text=True,
                capture_output=True
        )
        all_consts = glob(r'*.con') + glob(r'*.const')
        new_consts = [file for file in all_consts if file not in existing_cons]
        load_all = input(f'Load all constraint files [{', '.join(all_consts)}] (Y) or only new ones [{', '.join(new_consts)}]? (N) [Y]: ') or 'Y'
        if load_all.upper() == 'Y':
            const_files = all_consts
        else:
            const_files = new_consts
        print('Constraints file(s) loaded and will be added to masterfiles. Remember to check XDConstraints.out!')
        XDConOut = Path('XDConstraints.out')
        XDConOut.write_text(XDCon.stdout)

    elif Consts.upper() == 'W' and XDCon_path:
        print('Failed to locate and run XDConstraints.')
        print('Please ensure path to XDConstraints are correctly set (XDCon_path).')
        print('Alternatively run XDConstraints and add constraints manually.')
        print('Continuing with no constraint files.')
        const_files = []

    elif Consts.upper() == 'W':
        print('Please set the path to XDConstraints (XDCon_path).')
        print('Continuing with no constraint files.')
        const_files = []

    else:
        const_files = []

    # Checking loaded files
    print('\nFiles loaded:\nMaster:')
    print(''.join(oMas[:4]), end='')
    print('\nInst:\n', *inst_file, sep = '')
    print(f'{f'\nCons:\n[{', '.join(set(const_files))}]' if const_files else '\nNo constraints loaded.'}\n', sep = '')

    # Ensuring original masterfile is saved
    if not os.path.exists('xd00.mas'):
        os.rename(oMasName, 'xd00.mas')
    else:
        print(f'Original .mas ({oMasName}) file was not renamed/saved (xd00.mas already exists).')

    return oMas, inst_file, list(set(const_files))

class Atom():
    def __init__(self, atom, atom0, ax1, atom1, atom2, ax2, rl, tp, tbl, kap, lmx, sitesym = ' ', chemcon = ' '):
        self.atom = atom[1]
        self.coordsys = f'{atom0[1]:10}{ax1[1]:3}{atom1[1]:9}{atom2[1]:9}{ax2[1]:4}{rl[1]:4}'
        self.tp = tp[1]
        self.tbl = tbl[1]
        self.kap = kap[1]
        self.lmx = lmx[1]
        self.sitesym = sitesym[1]
        self.chemcon = '' if not chemcon else chemcon[1]
        self.instruction = []
        self.element = None

    def atom_table(self):
        chemcon = self.chemcon if 'CC' in self.instruction else ' '
        return f'{self.atom:9}{self.coordsys:39}{self.tp:3}{self.tbl:4}{self.kap:4}{self.lmx:3}{self.sitesym:6}{chemcon:5}'

    def key_table(self):
        xyz = '000'
        U = ' 000000 0000000000 000000000000000'
        MMs = {'M': '00',
               'D': '000',
               'Q': '00000',
               'O': '0000000',
               'H': '000000000'
               }

        # Handle hydrogen specific instructions correctly
        instructions = [ins if ins not in ['HXYZ', 'HU1'] else ins[1:] for ins in self.instruction]

        # Position
        if 'XYZ' in instructions: xyz = '111'

        # Thermal parameters
        if 'U1' in instructions: U = re.sub(' 0', ' 1', U, count = 1)
        if 'U2' in instructions: U = re.sub(' 000000 ', ' 111111 ', U, count = 1)
        if 'U3' in instructions: U = re.sub(' 0000000000 ', ' 1111111111 ', U)
        if 'U4' in instructions: U = re.sub(' 000000000000000', ' 111111111111111', U)
        # print(f'U: {U}')

        # Multipoles
        if 'M' in instructions: MMs['M'] = self.MPs['M']
        if 'D' in instructions: MMs['D'] = self.MPs['D']
        if 'Q' in instructions: MMs['Q'] = self.MPs['Q']
        if 'O' in instructions: MMs['O'] = self.MPs['O']
        if 'H' in instructions: MMs['H'] = self.MPs['H']

        Multipoles = ' '.join(str(mp) for mp in MMs.values())
        self.instruction = []
        return f'{self.atom:8}{xyz}{U} {Multipoles}'

    def set_pseudosymm(self, MPs):
        self.MPs = MPs

def clean_mas(masterfile):
    masterfile = deepcopy(masterfile)
    starters = ['ATOM     ATOM0    AX1', 'KEY']
    endings =  ['END ATOM', 'EXTCN ']
    new_master = []
    skip = False
    for line in masterfile:
        if re.match(r'!?include\s+', line, re.IGNORECASE):
            continue

        if re.match(r'SELECT\s+\*?model', line):
            new_master.append(re.sub(r'SELECT\s+\*?model(\s+\d){4}', 'SELECT  *model 0 0 0 0', line))
            continue

        if re.match(r'FOUR\s+fmod1', line):
            new_master.append(re.sub(r'FOUR\s+fmod1(\s+-?\d){4}\s+fmod2(\s+-?\d){4}',
                                      'FOUR     fmod1 0 0 0 0  fmod2 -1 0 0 0', line))
            continue
        
        if re.match(scale_regex, line):
            new_master.append(re.sub(r'[01]+', lambda m: '0' * len(m.group()), line))
            continue

        if any(line.startswith(marker) for marker in starters):
            new_master.append(line)
            skip = True
            continue

        if any(line.startswith(marker) for marker in endings):
            new_master.append(line)
            skip = False
            continue

        if not skip:
            new_master.append(line)

    with open('clean.mas', 'w') as cleanmas:
        cleanmas.writelines(new_master)

    return new_master

def get_atoms(masterfile):
    atomdict = {}
    scatdict = {}
    dummyatoms = []
    type_counter = 1
    for line in masterfile:

## Scat table (to get atom type)
        if scat_table_regex.search(line):
            scat = scat_table_regex.match(line).groupdict()
            if scat['SCAT'] in scatdict:
                scatdict[scat['SCAT']] += str(type_counter)
                type_counter += 1

            else:
                scatdict[scat['SCAT']] = str(type_counter)
                type_counter += 1

## Atom table
        if atom_table_regex.search(line):
            atom = atom_table_regex.match(line).groupdict()
            if atom['CHEMCON'] == None: atom['CHEMCON'] = ''
            atomdict[atom['ATOM']] = Atom(*atom.items())

        if dum_regex.search(line):
            dummyatoms.append(line)

## Key table (get MP pseudo symmetry)
        if key_table_regex.search(line):
            keys = key_table_regex.match(line).groupdict()
            MPs = {}
            for mp in ['M', 'D', 'Q', 'O', 'H']: MPs[mp] = keys.get(mp)
            atomdict[keys['ATOM']].set_pseudosymm(MPs)

        elif line == 'END KEY\n':
            break

## Finishing atoms (adding type)
    for atom in atomdict.values():
        for elemtype in scatdict.items():
            if atom.tbl in elemtype[1]:
                atom.element = elemtype[0]
                break

    return atomdict, dummyatoms

def get_instructions(line):
    dictionary = {}
    instruction = line.split()
    for i in instruction:
        re.sub(r'[\[\];]',' ', i).split()
        dictionary[re.sub(r'[\[\];]',' ', i).split()[0]] = re.sub(r'[\[\];]',' ', i).split()[1:]
    return dictionary

def get_other_values(masterfile:list):
    """Extract additional information from masterfile.

    Read masterfile to get random values such as:
    Number of scale factors
    
    Args:
        masterfile (list): Read-in masterfile

    Returns:
        scales (int): Number of scalefactors
        type_kappas (list): List of matching atom types and numbers, 
            and corresponding kappa indices
    """
    scales = 0
    type_kappas = []

    for line in masterfile:
        # Getting the number of scale factors
        if re.match(scale_regex, line):
            scales = len(re.search(r'[01]+', line).group())

        # Getting atom types + numbering
        if scat_table_regex.search(line):
            a_type = scat_table_regex.match(line).groupdict()['SCAT']
            type_kappas.append([a_type, len(type_kappas) + 1])

        # Adding kappa indexes to type_kappas
        if atom_table_regex.search(line):
            atom_line = atom_table_regex.match(line).groupdict()
            for scat_entry in type_kappas:
                if int(atom_line['TBL']) == scat_entry[1] and atom_line['KAP'] not in scat_entry[2:]:
                    scat_entry.append(atom_line['KAP'])

    return scales, type_kappas

def instruct_atoms(inst_line, atomdict, maxU = 0):

    other_ins = {}
    # Only scale refinement
    if all([inst in ['SCALE', 'CON', 'CC'] for inst in inst_line]):
        scales = {'SCALE': inst_line['SCALE'] or ['A']}
        if 'CON' in inst_line:
            scales['CON'] = inst_line['CON'] or []
        return True, scales, maxU
    
    # Handle non-atom specific entries
    for ins in [key for key in inst_line.keys() if key in ['SCALE', 'SIGOBS', 'SINTHL', 'KAPPA', 'KAPPAP', 'CON']]:
        if ins in ['SCALE', 'KAPPA', 'KAPPAP']:
            inst_line[ins] = inst_line[ins] or ['A']
        other_ins[ins] = inst_line.pop(ins)

    maxU_test = [U for U in inst_line.keys() if re.match(r'H?U\d', U)]
    maxU = int(max(maxU_test)[-1]) if (maxU_test and int(max(maxU_test)[-1]) > maxU) else maxU

    # Prepare element specific instructions
    inst_line = {
        ins: (
            'H' if ins in ('HXYZ', 'HU1')
            else v + ['!H'] if ins in ('XYZ', 'U2', 'U3', 'U4') and v
            else ['!H'] if ins in ('XYZ', 'U2', 'U3', 'U4')
            else ['A'] if not v
            else v
        )
        for ins, v in inst_line.items()
    }
    print(inst_line)

    # Handling element specific instructions
    for atom in atomdict.values():
        for ins in [key for key in inst_line.keys()]:
            inst = inst_line[ins]
            if '!' + atom.atom in inst or '!' + atom.element in inst: continue

            if (
                atom.atom in inst or atom.element in inst  # Element specified in instruction
                or all(s.startswith('!') for s in inst)    # All entries are negative (but not current element)
                or inst == ['A']                           # Instruction marked for all elements
            ):
                atom.instruction.append(ins)

    return False, other_ins, maxU

def write_kappa_entry(type_kappas, otherins):
    """
    Docstring for write_kappa_entry
    
    :param type_kappas: Description
    :param otherins: Description
    """

    # Collect lines by kappa index
    kappa_lines = defaultdict(list)

    def matches_selector(name, number, selector):
        if selector == ['A']:
            return True
        positives, negatives = set(), set()
        for s in selector:
            if s.startswith('!'):
                negatives.add(s[1:])
            else:
                positives.add(s)
        if positives and (name in positives or number in positives):
            return True
        if negatives and not (name in negatives or number in negatives):
            return True
        return False

    # Build lines
    for entry in type_kappas:
        name, number, *kappa_indices = entry
        bits = ['0']*6
        if 'KAPPA' in otherins and matches_selector(name, number, otherins['KAPPA']):
            bits[0] = '1'
        if 'KAPPAP' in otherins and matches_selector(name, number, otherins['KAPPAP']):
            bits[1] = '1'
        bitstring = ''.join(bits)

        for kidx in kappa_indices:
            kidx_int = int(kidx)
            kappa_lines[kidx_int].append((name, bitstring))

    # Write to file sorted by kappa index
    result = []
    for kidx in sorted(kappa_lines.keys()):
        for name, bitstring in kappa_lines[kidx]:
            result.append(f'KAPPA   {bitstring}\t! {name}\n')

    return result

def generate_con_includes(cons, otherins):
    """
    Generate INCLUDE / !INCLUDE lines for constraint files.
    Implements the 4 stated rules.
    """

    result = []
    # --- CASE 1: No cons at all ---
    if not cons:
        return result

    # Normalize cons into (base name → full file)
    con_map = {c.split('.')[0]: c for c in cons}
    

    # --- CASE 2: CON not present → no con files included ---
    if 'CON' not in otherins:
        return [f'!INCLUDE {c}' for c in cons]

    # Now CON is present
    con_values = otherins['CON']

    # If CON present but empty list -> include everything
    if not con_values:
        return [f'INCLUDE {c}' for c in cons]

    # Collect explicit includes/excludes
    includes = set()
    excludes = set()

    for entry in con_values:
        stripped = entry.lstrip('!')
        base = stripped.split('.')[0]

        if entry.startswith('!'):
            excludes.add(base)
        else:
            includes.add(base)

    # --- CASE 3: Positive inclusion (only files in includes) ---
    if includes:
        # Ignore any excludes for this rule — positive list overrides
        for base, full in con_map.items():
            if base in includes:
                result.append(f'INCLUDE {full}')
            else:
                result.append(f'!INCLUDE {full}')

        return sorted(result, key = lambda s: s.startswith('!'))

    # --- CASE 4: Negative inclusion only (no positive includes) ---
    # Include everything *except* excluded ones
    for base, full in con_map.items():
        if base in excludes:
            result.append(f'!INCLUDE {full}')
        else:
            result.append(f'INCLUDE {full}')

    return sorted(result, key = lambda s: s.startswith('!'))

def write_masterfile(masterfile:list, atomdict:dict, dummyatoms:list, cons:list, number:int, maxU:int, only_scale:bool = False, otherins:dict = {}):
    """
    Docstring for write_masterfile
    
    Args:
        masterfile (list): Input masterfile
        atomdict (dict): Dictionary of all atoms in masterfile
        dummyatoms (list): List of dummy atoms (read from masterfile)
        cons (list): List of constraint files from input
        number (int): Current instruction number
        maxU (int): Max U-parameter so far
        only_scale (bool): Only refining scalefactors
        otherins (dict): Dictionary of non-atom specific instructions
    """
    ## Prepping
    number = f'{number:02}'
    scales, type_kappas = get_other_values(masterfile)
    cleanmas = clean_mas(masterfile)

    ## Writing
    with open(f'xd{number}.mas', 'w') as new_master:
        for i, line in enumerate(cleanmas):
            # End of XDLSM -> Leave the rest of xd.mas untouched
            if line == '   END XDLSM\n':
                new_master.writelines(cleanmas[i:])
                break

            # Cycles (always 30)
            elif re.match(cycles_regex, line):
                if only_scale:
                    line = re.sub(r'cycle\s+-?\d+\s+', 'cycle -30 ', line)
                else:
                    line = re.sub(r'cycle\s+-?\d+\s+', 'cycle  30 ', line)

            # Model (always uses lmx = 4)
            elif re.match(r'SELECT\s+\*?model\s+-?\d\s+-?\d\s+\d\s+\d', line):
                line = re.sub(r'\*?model\s+-?\d\s+-?\d\s+\d\s+\d', f'*model 4 {maxU} 1 0', line)

            # FOUR line (always uses lmx = 4/-1 (fmod1/fmod2)) and constraints if present
            elif re.match(four_regex, line):
                line = f'FOUR     fmod1 4 {maxU} 0 0  fmod2 -1 {maxU} 0 0\n'
                new_master.write(line)

                constraints = generate_con_includes(cons, otherins)
                for con in constraints:
                    new_master.write(con + '\n')
                continue

            # Scale factors
            elif re.match(scale_regex, line) and 'SCALE' in otherins:
                if 'A' in otherins['SCALE']:
                    line = re.sub(r'[01]+', '1' * scales, line)
                    new_master.write(line)
                    continue
                
                if any('!' in scale for scale in otherins['SCALE']):
                    indexes = list(range(1, scales + 1))
                    for scale in otherins['SCALE']:
                        indexes.remove(int(scale[1:]))
                else:
                    indexes = []
                    for scale in otherins['SCALE']:
                        indexes.append(int(scale))

                sfs = '0' * scales
                for i in indexes: sfs = sfs[:i-1] + '1' + sfs[i:]
                line = re.sub(r'[01]+', sfs, line)
                new_master.write(line)
                continue

            # SINTHL and SIGOBS
            elif re.match(skip_regex, line):
                if 'SINTHL' in otherins:
                    lower, upper = otherins['SINTHL']
                    line = re.sub(r'[\*]{0,1}sinthl\s+[d\d.]+\s+[d\d.]+', f'*sinthl {lower} {upper}', line)
                
                if 'SIGOBS' in otherins:
                    limit = float(otherins['SIGOBS'][0])
                    line = re.sub(r'[\*]{0,1}sigobs\s+[d\d.]+', f'*sigobs {limit:.1f}', line)

            # Atom table
            elif line.startswith('ATOM     ATOM0    AX1'):
                new_master.write(line)
                [new_master.write(atom.atom_table() + '\n') for atom in atomdict.values()]
                [new_master.write(dum) for dum in dummyatoms]
                continue

            # Key table (and kappas)
            elif line.startswith('KEY'):
                new_master.write(line)
                [new_master.write(atom.key_table() + '\n') for atom in atomdict.values()]
                kappa_list = write_kappa_entry(type_kappas, otherins)
                [new_master.write(kappa_entry) for kappa_entry in kappa_list]

                continue

            new_master.write(line)

def main():
    """
    Main function :))
    """
    # Getting files
    mas, inst, cons = get_files()

    # Defining constants
    maxU = 2  # Maximum TP value (so far), starts at anisotropic (2).

    # Gathering information
    atoms, dummys = get_atoms(mas)

    for number, line in enumerate(inst, start = 1):
        insts = get_instructions(line)
        onlyscale, other_ins, maxU = instruct_atoms(insts, atoms, maxU)
        write_masterfile(mas, atoms, dummys, cons, number, maxU, onlyscale, other_ins)

    # Cleaning up
    os.remove(r'clean.mas')

if __name__ == '__main__':
    main()
