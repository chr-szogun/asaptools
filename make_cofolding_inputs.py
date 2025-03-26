import json
from rdkit import Chem
from pymol import cmd
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument("-n", "--name", type=str, help="name of the protein", default="protein")
parser.add_argument("-p", "--proteinpdb", type=str, help="pdb file containing the protein", default="protein.pdb")
parser.add_argument("-l", "--ligandsdf", type=str, help="sdf containing the ligands", default="ligands.sdf")
args = parser.parse_args()

def read_ligands_sdf(path)->Chem.SDMolSupplier:
    """reads the ligand.sdf and returns the rdkit SDMolSupplier object"""
    mols = Chem.SDMolSupplier(path)
    return mols

def read_proteinseq(path)->list:
    """reads the protein.pdb and returns the fasta sequence as a string"""
    cmd.delete("all")
    cmd.load(path,"my_prot")
    cmd.select("polymer")
    seq_raw = cmd.get_fastastr("sel").split(">my_prot")
    seq = [re.sub("_[A-Z]\\n","",l) for l in seq_raw]
    seq = [s.replace("\n","") for s in seq if not s=='']
    print(f"found {len(seq)} protein sequence(s) of length(s) {[len(s) for s in seq]} in {path}.")
    return seq


def write_boltz_input(proteinseq,pname,ligsmiles,ligname):
    """takes the protein sequence, protein name, ligand SMILES, and ligand name
    and writes a boltz fasta input from scratch"""

    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXY"
    with open(f"boltz_{pname}_{ligname}.fasta", "w+") as o:
        for i in range(len(proteinseq)):
            o.write(f">{alphabet[i]}|protein\n{proteinseq[i]}\n")
        o.write(f">Z|smiles\n{ligsmiles}")


def write_chai_input(proteinseq,pname,ligsmiles,ligname):
    """takes the protein sequence, protein name, ligand SMILES, and ligand name
    and writes a chai fasta input using the template"""

    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXY"
    with open(f"chai_{pname}_{ligname}.fasta", "w+") as o:
        for i in range(len(proteinseq)):
            o.write(f">protein|name={pname}_{alphabet[i]}\n{proteinseq[i]}\n")
        o.write(f">smiles|name={ligname}\n{ligsmiles}")


def write_af_input(proteinseq,pname,ligsmiles,ligname):
    """takes the protein sequence, protein name, ligand SMILES, and ligand name
    and writes an alphafold3 json input using the template"""

    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXY"
    d = {"name":f"{pname}_{ligname} cofold run","modelSeeds":[42],
         "sequences":[],"dialect":"alphafold3","version":3}
    
    for i in range(len(proteinseq)):
        proteinblock = {"protein":{"id":alphabet[i], "sequence":proteinseq[i]}}
        d["sequences"].append(proteinblock)
    d["sequences"].append({"ligand":{"id":"Z","smiles":ligsmiles}})

    with open(f"af_{pname}_{ligname}.json", "w+") as o:
        json.dump(d,o,indent=4)


def main(pname,proteinpdb,ligandsdf):
    mols = read_ligands_sdf(ligandsdf)
    proteinseq = read_proteinseq(proteinpdb)
    for mol in mols:
        ligname = mol.GetProp("_Name").replace(" ","")
        ligsmiles = Chem.MolToSmiles(mol)
        write_boltz_input(proteinseq,pname,ligsmiles,ligname)
        write_chai_input(proteinseq,pname,ligsmiles,ligname)
        write_af_input(proteinseq,pname,ligsmiles,ligname)

main(args.name, args.proteinpdb, args.ligandsdf)