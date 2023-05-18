import os

for Li in range(5):
    for Bi in ['B0', 'B1', 'B2']:
    # for Bi in ['Bnot_hydro']:

        command_str = f"mv mix_L{Li}_Ma0_{Bi} mix_L{Li}_Ma10_{Bi}"
        os.system(command_str)

        print(f"{command_str}")

