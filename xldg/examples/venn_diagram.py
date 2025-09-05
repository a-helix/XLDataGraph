import os
from xldg.data import Path, MeroX
from xldg.graphics import VennConfig, Venn2, Venn3

if __name__ == "__main__":
    # Set up working directory
    cwd = os.path.join(os.getcwd(), 'examples', 'files')
    print(f"Current Working Directory: {cwd}\n")

    # 1: Load crosslink data files
    crosslink_files = Path.list_given_type_files(cwd, 'zhrm')
    crosslinks = MeroX.load_data(crosslink_files, 'DSBU')
    print("#1: Crosslink data loaded\n")

    # 2: Configure Venn diagram settings
    config = VennConfig('Batch1', 'Batch2', 'Batch3', title='Title')
    print("#2: Venn diagram configuration created\n")

    # 3: Build and save Venn2 diagram
    venn2 = Venn2(crosslinks[0], crosslinks[1], config)
    venn2_path = os.path.join(cwd, 'venn2.svg')
    venn2.save(venn2_path)
    print(f"#3: Venn2 diagram saved to {venn2_path}\n")

    # 4: Build and save Venn3 diagram
    venn3 = Venn3(crosslinks[0], crosslinks[1], crosslinks[2], config)
    venn3_path = os.path.join(cwd, 'venn3.svg')
    venn3.save(venn3_path)
    print(f"#4: Venn3 diagram saved to {venn3_path}\n")
