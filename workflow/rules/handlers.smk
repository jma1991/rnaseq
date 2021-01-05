# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

class colors:
    HEADER = '\033[95m'
    PASS = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

# Print workflow information upon starting
onstart:
    print(colors.BOLD + colors.OKBLUE + "# Design" + colors.ENDC)
    print(colors.OKBLUE + "# Conditions:", ", ".join(CONDITIONS) + colors.ENDC)
    print(colors.OKBLUE + "# Contrasts:", ", ".join(CONTRASTS) + colors.ENDC)
    print(colors.OKBLUE + "# Samples:", ", ".join(SAMPLES) + colors.ENDC)
    print(colors.OKBLUE + "# Runs:", ", ".join(RUNS) + colors.ENDC)
    print(colors.OKBLUE + "# Genome:", GENOME + colors.ENDC)

# Print message when workg success workflow
onsuccess:
    print(colors.BOLD + colors.PASS + "Workflow finished!" + colors.ENDC)

# Print message when workflow encounters an error
onerror:
    print(colors.BOLD + colors.FAIL + "An error occurred!" + colors.ENDC)