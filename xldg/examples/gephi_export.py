import os
from xldg.data import Path, MeroX, CrossLink, ProteinChain

if __name__ == "__main__":
    # Set up working directory
    cwd = os.path.join(os.getcwd(), 'examples', 'files')
    print(f"Current Working Directory: {cwd}\n")

    # 1: Load Protein Chain Data
    pcd_path = os.path.join(cwd, "two_monomers_complex.pcd")
    pcd = ProteinChain.load_data(pcd_path)
    print("#1: Protein chain data loaded\n")

    # 2: Load and combine crosslink (zhrm) data
    zhrm_folder_path = Path.list_given_type_files(cwd, 'zhrm')
    linker = 'DSBU'
    folder_content = MeroX.load_data(zhrm_folder_path, linker)
    combined_dataset = CrossLink.combine_all(folder_content)
    print("#2: Crosslink data loaded and combined\n")

    # 3: Export crosslinks for Gephi visualization
    save_folder = os.path.join(cwd, 'results')

    # 3.1: Export protein-protein interactions (PPIs)
    ppi_file = os.path.join(save_folder, "PPIs.gexf")
    combined_dataset.export_ppis_for_gephi(save_folder, "PPIs.gexf", pcd)
    print(f"\t#3.1: Protein-protein interactions exported to {ppi_file}")

    # 3.2: Export amino acid interactions (AAIs)
    aai_file = os.path.join(save_folder, "AAIs.gexf")
    combined_dataset.export_aais_for_gephi(save_folder, "AAIs.gexf", pcd)
    print(f"\t#3.2: Amino acid interactions exported to {aai_file}")

    print("#3: Done\n")
