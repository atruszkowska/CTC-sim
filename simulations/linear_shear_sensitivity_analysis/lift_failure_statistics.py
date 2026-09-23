from pathlib import Path
import re
import pandas as pd


def process_bond_failures(directory=".", output_file="single_cell_lift_coefficient_bond_failure_statistics.csv"):
    directory_path = Path(directory)
    results = []

    # Matches strictly:
    # 1. bond_failures<index>_base.txt
    # 2. bond_failures<index>_Single_cell_lift_coefficient_<value>.txt
    pattern = re.compile(
        r"^bond_failures(?P<index>\d+)(?:_(?P<is_base>base)|_Single_cell_lift_coefficient_(?P<coefficient>[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?))\.txt$"
    )

    for file_path in directory_path.glob("bond_failures*.txt"):
        match = pattern.match(file_path.name)
        if not match:
            continue

        file_info = match.groupdict()

        # Handle lift coefficient and base case
        if file_info["is_base"]:
            coefficient = None
            property_name = "base"
        else:
            coefficient = float(file_info["coefficient"])
            property_name = "Single_cell_lift_coefficient"

        # Read the file data
        try:
            df = pd.read_csv(file_path, header=None)
            values = df.iloc[:, 0]
            total_count = len(values)

            if total_count == 0:
                continue

            # Compute percentages for 0, 1, and 2
            counts = values.value_counts()
            pct_0 = (counts.get(0, 0) / total_count) * 100
            pct_1 = (counts.get(1, 0) / total_count) * 100
            pct_2 = (counts.get(2, 0) / total_count) * 100

            results.append(
                {
                    "filename": file_path.name,
                    "index": int(file_info["index"]),
                    "property": property_name,
                    "value": (
                        "base"
                        if property_name == "base"
                        else str(coefficient)
                    ),
                    "pct_0": pct_0,
                    "pct_1": pct_1,
                    "pct_2": pct_2,
                    "total_bonds": total_count,
                }
            )

        except Exception as e:
            print(f"Error reading {file_path.name}: {e}")

    # Convert results to DataFrame and save to CSV
    if results:
        summary_df = pd.DataFrame(results)
        summary_df.sort_values(
            by=["property", "index"],
            inplace=True,
        )

        summary_df.to_csv(output_file, index=False)

        print(
            f"Successfully processed {len(results)} files. "
            f"Summary saved to '{output_file}'."
        )
    else:
        print("No matching bond_failures files found.")


if __name__ == "__main__":
    process_bond_failures()
