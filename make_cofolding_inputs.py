import yaml
from rdkit import Chem
from pymol import cmd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-n", "--name", type=str, help="name of the protein", default="protein")
parser.add_argument("-p", "--proteinpdb", type=str, help="pdb file containing the protein", default="protein.pdb")
parser.add_argument("-l", "--ligandsdf", type=str, help="sdf containing the ligands", default="ligands.sdf")
args = parser.parse_args()

#def read_ligands_yml(path):
#    with open(path, "r") as f:
#        ligyml = yaml.safe_load(f)
#    return ligyml

def read_ligands_sdf(path)->Chem.SDMolSupplier:
    """reads the ligand.sdf and returns the rdkit SDMolSupplier object"""
    mols = Chem.SDMolSupplier(path)
    return mols

def read_proteinseq(path)->str:
    """reads the protein.pdb and returns the fasta sequence as a string"""
    cmd.load(path,"my_prot")
    cmd.select("polymer")
    seq_raw = cmd.get_fastastr("sel").split()[1:]
    seq = "".join(seq_raw)
    return seq

def write_boltz_input(proteinseq,pname,ligsmiles,ligname):
    """takes the protein sequence, protein name, ligand SMILES, and ligand name
    and writes a boltz fasta input using the template"""
    with open("templates/boltz_template.fasta", "r") as f:
        with open(f"boltz_{pname}_{ligname}.fasta", "w+") as o:
            for line in f:
                o.write(line.replace("PROTEINFASTA",proteinseq).replace("LIGANDSMILES",ligsmiles))


def write_chai_input(proteinseq,pname,ligsmiles,ligname):
    """takes the protein sequence, protein name, ligand SMILES, and ligand name
    and writes a chai fasta input using the template"""
    with open("templates/chai_template.fasta", "r") as f:
        with open(f"chai_{pname}_{ligname}.fasta", "w+") as o:
            for line in f:
                o.write(line.replace("PROTEINFASTA",proteinseq).replace("LIGANDSMILES",ligsmiles))


def write_af_input(proteinseq,pname,ligsmiles,ligname):
    """takes the protein sequence, protein name, ligand SMILES, and ligand name
    and writes an alphafold3 json input using the template"""
    with open("templates/alphafold_template.json", "r") as f:
        with open(f"af_{pname}_{ligname}.json", "w+") as o:
            for line in f:
                o.write(line.replace("PROTEINFASTA",proteinseq).replace("LIGANDSMILES",ligsmiles))


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