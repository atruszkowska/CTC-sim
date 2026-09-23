import math

results_file = "study_results_parametric.txt"
target_velocity = 2.5e-6 
top_n = 10

entries = []
lines = [line.rstrip("\n") for line in open(results_file)]

i = 0
while i < len(lines):
    if lines[i].startswith("Parameters:"):
        # Extract the 4-line block
        block = lines[i:i+4]
        params = block[0].split("Parameters:")[1].strip()

        # Extract velocity from line 4 in the block
        velocity_str = block[3].split(":")[1].strip()
        velocity = float(velocity_str)

        diff = abs(velocity - target_velocity)/target_velocity*100
        entries.append((diff, params, velocity, block))

        i += 4  # move to next block
    else:
        i += 1

# Sort by closeness to target velocity
entries.sort(key=lambda x: x[0])

# Take top N
top_entries = entries[:top_n]

print(f"Top {top_n} entries closest to target velocity {target_velocity}:\n")
for rank, (diff, params, velocity, block) in enumerate(top_entries, start=1):
    print(f"--- Entry #{rank} ---")
    for bline in block:
        print(bline)
    print(f"Velocity relative error, %: {diff:.3e}")
    print()

