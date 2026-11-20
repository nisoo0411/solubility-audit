"""Regenerate six descriptors and add RDKit scaffold groups for the benchmark."""

import argparse
import csv
from io import StringIO
from pathlib import Path
from urllib.request import urlopen

from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

DATA_URL = (
    "https://raw.githubusercontent.com/deepchem/deepchem/master/"
    "datasets/delaney-processed.csv"
)
FEATURES = {
    "Minimum Degree": lambda mol: min(atom.GetDegree() for atom in mol.GetAtoms()),
    "Molecular Weight": Descriptors.MolWt,
    "Number of H-Bond Donors": Lipinski.NumHDonors,
    "Number of Rings": Lipinski.RingCount,
    "Number of Rotatable Bonds": Lipinski.NumRotatableBonds,
    "Polar Surface Area": rdMolDescriptors.CalcTPSA,
}


def featurize(rows):
    for row in rows:
        mol = Chem.MolFromSmiles(row["smiles"])
        if not mol or not mol.GetNumAtoms():
            raise ValueError(f"Invalid SMILES for {row['Compound ID']!r}")
        scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
        yield {
            "Compound ID": row["Compound ID"],
            "smiles": Chem.MolToSmiles(mol),
            # Group all acyclic structures together so the split also tests transfer to them.
            "scaffold": scaffold or "ACYCLIC",
            "measured log solubility in mols per litre": row[
                "measured log solubility in mols per litre"
            ],
            **{name: calculate(mol) for name, calculate in FEATURES.items()},
        }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, help="Local Delaney CSV")
    parser.add_argument("--output", type=Path, default=Path("rdkit-esol.csv"))
    args = parser.parse_args()
    if args.input:
        raw = args.input.read_text(encoding="utf-8")
    else:
        with urlopen(DATA_URL, timeout=30) as response:
            raw = response.read().decode("utf-8")
    rows = list(featurize(csv.DictReader(StringIO(raw))))
    with args.output.open("w", encoding="utf-8", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
