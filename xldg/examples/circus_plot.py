import os
from xldg.data import Path, MeroX, Domain, Fasta, CrossLink
from xldg.graphics import CircosConfig, Circos

if __name__ == "__main__":
    # Set up working directory
    cwd = os.path.join(os.getcwd(), 'examples', 'files')
    print(f"Current Working Directory: {cwd}\n")

    # 1: Load FASTA dataset
    fasta_path = os.path.join(cwd, 'example.fasta')
    fasta_format = 'Custom'
    fasta_dataset = Fasta.load_data(fasta_path, fasta_format)
    print("#1: FASTA data loaded\n")

    # 2: Load domain files
    domain_path = Path.list_given_type_files(cwd, 'dmn')
    domains = Domain.load_data(domain_path)
    print("#2: Domain data loaded\n")

    # 3: Load crosslink (zhrm) data
    zhrm_folder_path = Path.list_given_type_files(cwd, 'zhrm')
    linker = 'DSBU'
    folder_content = MeroX.load_data(zhrm_folder_path, linker)
    print("#3: Crosslink (zhrm) data loaded\n")

    # 4: Combine and filter crosslink data
    combined_data = CrossLink.combine_all(folder_content)
    filtered_data = CrossLink.filter_by_score(combined_data, min_score=30)
    print("#4: Crosslink data combined and filtered\n")

    # 5: Generate Circos plot
    save_path = os.path.join(cwd, 'results', "circos_basic.svg")
    config = CircosConfig(fasta_dataset, domains)
    circos = Circos(combined_data, config)
    circos.save(save_path)
    print(f"#5: Circos plot saved to {save_path}\n")
