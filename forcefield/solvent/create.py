from rdkit import Chem
from rdkit.Chem import AllChem

# Read SMILES from file manually (since we see it clearly)
mol = Chem.MolFromSmiles("CCCCCCCCCCCCCCCC")
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol)
AllChem.MMFFOptimizeMolecule(mol)

# Set residue info for each atom
for i, atom in enumerate(mol.GetAtoms(), start=1):
    info = Chem.AtomPDBResidueInfo()
    info.SetName(f"{atom.GetSymbol():>2s}")  # atom name
    info.SetResidueName("HEX")               # residue name
    info.SetResidueNumber(1)                 # residue number
    info.SetChainId("A")                     # optional
    atom.SetMonomerInfo(info)

Chem.MolToPDBFile(mol, "hexadecane-raw.pdb")
