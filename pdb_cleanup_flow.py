from pymol import cmd
import glob
from rdkit import Chem
from rdkit.Chem import AllChem
from asapdiscovery.data.schema.complex import Complex, PreppedComplex
import os
import json
import re

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-r", "--referenceligand", type=str,
                     help="name of the reference ligand that forms the complex all others should be aligned to")
parser.add_argument("-p", "--pathtopdbs", type=str,
                     help="path to the pdb files")
args = parser.parse_args()

class complex_aligner():

    def __init__(self,ref):
        self.ref = ref

        self.align_all_complexes(self.ref)

    def align_to_ref(self,to_align:str, ref:str):
        cmd.delete("all")
        name = to_align.replace("_model_0.pdb","")
        cmd.load(to_align, name)
        cmd.load(ref, "ref")
        cmd.align(name, "ref")
        cmd.save(f"{name}_aligned.pdb",name)

    def align_all_complexes(self,ref:str):
        for pdb in glob.glob("*_model_0.pdb"):
            ref_pdb = "_".join(pdb.split()[:1][2:])
            self.align_to_ref(pdb, ref_pdb)


class LigBondOrder_fixer():

    def __init__(self):
        self.main()


    def read_ligand_from_PDB(self,pdb_file)->Chem.rdchem.Mol:
        ligand_pdbblock = ""
        with open(pdb_file, "r") as pdb:
            for line in pdb:
                if line.startswith("HETATM"):
                    ligand_pdbblock += line
                elif line.startswith("CONECT"):
                    ligand_pdbblock += line
        ligand_3D = Chem.MolFromPDBBlock(ligand_pdbblock)
        #smi = Chem.MolToSmiles(ligand_3D)
        #ligand_2D = Chem.MolFromSmiles(smi)
        return ligand_3D

    def read_ligand_from_fasta(self,fasta_file)->Chem.rdchem.Mol:
        with open(fasta_file, "r") as fasta:
            smiles = fasta.readlines()[-1].strip()
        ligand = Chem.MolFromSmiles(smiles)
        return ligand

    def read_ligand_from_json(self,json_file)->Chem.rdchem.Mol:
        import json
        with open(json_file, "r") as f:
            d = json.load(f)
            smiles = d["sequences"][1]["ligand"]["smiles"]
        ligand = Chem.MolFromSmiles(smiles)
        return ligand

    def are_ligand_smiles_identical(self,ligand1:Chem.rdchem.Mol, ligand2:Chem.rdchem.Mol)->bool:
        smiles1 = Chem.CanonSmiles(Chem.MolToSmiles(ligand1))
        smiles2 = Chem.CanonSmiles(Chem.MolToSmiles(ligand2))
        return True if smiles1==smiles2 else False


    def assign_bond_order(self,ligand_to_correct:Chem.rdchem.Mol, template:Chem.rdchem.Mol)->Chem.rdchem.Mol:
        ligand_corrected = AllChem.AssignBondOrdersFromTemplate(template, ligand_to_correct)
        return ligand_corrected

    def rewrite_PDB(self,pdb_file_old, ligand, pdb_file_new):
        with open(pdb_file_new, "w+") as pdb_new:
            with open(pdb_file_old, "r") as pdb_old:
                for line in pdb_old:
                    if line.startswith("HETATM"):
                        break
                    pdb_new.write(line)
                pdbblock_lig = Chem.MolToPDBBlock(ligand)
                pdb_new.write(pdbblock_lig)

    def read_ligand_ref(self,path_in)->Chem.rdchem.Mol:
        if path_in.endswith("fasta"):
            ligand = self.read_ligand_from_fasta(path_in)
        elif path_in.endswith("json"):
            ligand = self.read_ligand_from_json(path_in)
        return ligand

    def run_single(self,pdb_in, ref_in, pdb_out):
        lig_pdb3D = self.read_ligand_from_PDB(pdb_in)
        lig_ref = self.read_ligand_ref(ref_in)
        if not self.are_ligand_smiles_identical(lig_pdb3D, lig_ref):
            lig_corr = self.assign_bond_order(lig_pdb3D, lig_ref)
            self.rewrite_PDB(pdb_in, lig_corr, pdb_out)
        else:
            self.rewrite_PDB(pdb_in, lig_pdb3D, pdb_out)

    def main(self):
        for pdb in glob.glob(f"*aligned.pdb"):
            ref_in = pdb.replace("_aligned.pdb",".fasta")
            outname = pdb.replace(".pdb", "_fixedlig.pdb")
            self.run_single(pdb, ref_in, outname)


class protein_prep():

    def __init__(self, ref):
        self.ref = ref

        self.run_proteinprep(self.ref)

    def run_proteinprep(self, ref):
        for pdb in glob.glob("*aligned_fixedlig.pdb"):
            ligname = pdb.split("_")[-3]
            protname = pdb.split("_")[-4]

            cplx = PreppedComplex.from_complex(Complex.from_pdb(pdb, target_kwargs={"target_name":protname},ligand_kwargs={"compound_name":ligname}))
            cplx.to_json_file(f"{protname}_{ligname}.json")

            if ligname == ref: #most potent
                cplx.target.to_pdb_file(f"{protname}_boltz_raw.pdb") 


class SDF_writer():

    def __init__(self):
        self.json_to_sdf()

    def json_to_sdf():
        sdfs = []

        for jsonfile in glob.glob("*.json"):

            print(jsonfile)
            name = jsonfile.split("_")[1].split(".json")[0]
            with open(jsonfile,"r") as f:
                sdf = json.load(f)["ligand"]["data"]
                sdf = re.sub("LIG\(.*\)",name, sdf)  # I dont know why, but asapdiscovery keeps replacing the ligname with stuff like this
                #sdf = sdf.replace("LIG(B-289)",name)
                sdfs.append(sdf)

        with open("allligs.sdf", "w+") as s:
            [s.write(sdf) for sdf in sdfs]


os.chdir(args.pathtopdbs)
complex_aligner()
LigBondOrder_fixer()
protein_prep()
SDF_writer()