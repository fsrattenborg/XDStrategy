# XDStrategy
_Original by Lennard Krause_
Generate successive xd.mas files based on instructions file.
To be used with XDRefine (by Lennard Krause) or alike for sequential refinements in XD.

## Setup
- Set path of XDConstraints (if to be used).
- Initial configuration of masterfile (see below).
- Run script in folder with input files.

## Input files
- XD masterfile.
- Instructions file (standard sequence writte by the script).
- Constraint files (optional).
    - Will be included using INCLUDE commands.

### Initial masterfile setup
- Configure local coordinate systems for all atoms.
- Set/check symmetry restrictions on coordinates, ADPs and multipoles (including pseudosymmetries).
- Set values (TP/KAP/LMX/etc.) in ATOM table. 
