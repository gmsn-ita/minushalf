import os
from minushalf.utils import ElectronicDistribution

for element in ElectronicDistribution:
    symbol = str(element)
    path = f'/home/henrique/Documents/ic_physics/minushalf/tests/fixtures/{symbol}/'
    with open(os.path.join(path, "INP"), "r") as inp:
        lines = inp.readlines()
        lines[0] = lines[0].rstrip('\n') + "#Comment 1\n"
        lines[1] = lines[1].rstrip('\n') + "#Comment2\n"
        lines[3] = lines[3].rstrip('\n') + "#Comment4\n"
        lines[2] = lines[2].rstrip('\n') + "#CommentComment3\n"

        new_lines = [f"#{symbol}\n", "#all electrons\n", *lines]
        with open(os.path.join(path, "INP_COMMENTED"), "w") as file:
            file.writelines(new_lines)
