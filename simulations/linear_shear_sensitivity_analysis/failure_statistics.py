from pathlib import Path
import re
import pandas as pd


def process_bond_failures(directory=".", output_file="bond_failure_statistics.csv"):
    directory_path = Path(directory)
    results = []

    # Regex to capture the prefix index, property name (optional), and value (optional)
    # Matches both standard files and '_base.txt' files
    pattern = re.compile(
        r"^bond_failures(?P<index>\d+)(?:_(?P<property>.*?))?_(?P<value>-?\d+\.?\d*(?:e[+-]?\d+)?|base)\.txt$"
    )

    for file_path in directory_path.glob("bond_failures*.txt"):
        match = pattern.match(file_path.name)
        if not match:
            continue

        file_info = match.groupdict()

        # Read the file data (assumes a single column or whitespace-separated values)
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
                    "property": file_info["property"]
                    if file_info["property"]
                    else "base",
                    "value": file_info["value"],
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
        summary_df.sort_values(by=["index", "property"], inplace=True)
        summary_df.to_csv(output_file, index=False)
        print(
            f"Successfully processed {len(results)} files. Summary saved to '{output_file}'."
        )
    else:
        print("No matching bond_failures files found.")


if __name__ == "__main__":
    process_bond_failures()
