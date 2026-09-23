from pathlib import Path
import sys

import pandas as pd
sys.path.insert(0, '../../scripts/')
import utils as ut
import simulation_utils as smut
from colors import *

# (1) Read Positions Test

# Build a small Dataframe
expected = pd.DataFrame(
        {
            "x":  [0.10, 1.10],
            "y":  [0.20, 1.20],
            "it": [0,    1],
            "ir": [0,    0],
        }
)

# Write it out as txt file
name = "Read positions file"
sample_path = Path("test_pos.txt")
expected.to_csv(sample_path, sep=" ", index=False, header=False)

# Read the generated temp file
pos_data = smut.read_pos(sample_path)

# Compare
ok = pos_data.equals(expected)
ut.test_pass(ok, name)

# (2) Read Bonds Test
records = pd.DataFrame([
        ("1-2", 10, 0),
        ("3-4", 10, 0),
        ("1-3", 20, 0),
        ("2-4", 20, 0),
    ], columns=["ij", "it", "ir"], 
)

expected = pd.DataFrame(
        {
            "ID-1":  ["1", "3", "1", "2"],
            "ID-2":  ["2", "4", "3", "4"],
            "it": [10, 10, 20, 20],
            "ir": [0, 0, 0, 0],
        }
        )

# Write out bond_data to file
name = "Read bonds file"
sample_path = Path("test_bond.txt")
records.to_csv(sample_path, sep=" ", index=False, header=False)

# Read the generated file
bond_data = smut.read_bonds("test_bond.txt")

# Compare
ok = bond_data.equals(expected)
ut.test_pass(ok, name)
