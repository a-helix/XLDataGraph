import os, time
from xldg.data import Path, MeroX, CrossLink, ProteinStructure, ProteinChain

if __name__ == "__main__":
    # Set up working directory
    cwd = os.path.join(os.getcwd(), 'examples', 'files')
    print(f"Current Working Directory: {cwd}\n")

    # 1: Load Protein Chain Data
    pcd_path = os.path.join(cwd, "two_monomers_complex.pcd")
    pcd = ProteinChain.load_data(pcd_path)
    print("#1: Protein chain data loaded\n")

    # 2: Load Protein Structure Data
    structure_path = os.path.join(cwd, "two_monomers_complex.cif")
    structure = ProteinStructure.load_data(structure_path)
    print("#2: Protein structure data loaded\n")

    # 3: Load and process crosslink data
    crosslink_files = Path.list_given_type_files(cwd, 'zhrm')
    linker = 'DSBU'
    crosslinks = MeroX.load_data(crosslink_files, linker)

    print("\tCombining all replicates and blanking crosslink counters for simplicity")
    crosslinks = CrossLink.combine_all(crosslinks)
    crosslinks = CrossLink.blank_replica(crosslinks)
    print("#3: Crosslink data loaded and processed\n")

    # 4: Export crosslinks to ChimeraX
    print("#4: Exporting crosslinks to ChimeraX")

    # 4.1: Export without protein structure
    experiment_name = 'IMP_crosslinks_without_structure'
    export_folder = os.path.join(cwd, 'results', experiment_name)
    print(f"\t#4.1: Exporting without protein structure → {export_folder}")
    crosslinks.export_for_chimerax(pcd, export_folder, experiment_name)

    # 4.2: Export with protein structure and distance constraints
    experiment_name = 'IMP_crosslinks_with_structure'
    export_folder = os.path.join(cwd, 'results', experiment_name)
    min_distance, max_distance = 2.0, 35.0
    print(f"\t#4.2: Exporting with protein structure and distance constraints → {export_folder}")
    crosslinks.export_for_chimerax(pcd, export_folder, experiment_name,
                                   protein_structure=structure,
                                   min_distance=min_distance,
                                   max_distance=max_distance)
    print("#4: Done\n")

    # 5: Predict crosslinks
    print("#5: Predicting crosslinks")

    # 5.1: Prediction using direct distances between CA atoms
    from_residues, to_residues = '{K', '{K'
    experiment_name = 'IMP_crosslinks_prediction_direct'
    export_folder = os.path.join(cwd, 'results', experiment_name)
    print(f"\t#5.1: Predicting crosslinks (direct CA distances) → {export_folder}")

    prediction = structure.predict_crosslinks(
        pcd,
        from_residues,
        to_residues,
        min_distance,
        max_distance,
        'pseudoDSBU'
    )
    prediction.export_for_chimerax(pcd, export_folder, experiment_name)

    # 5.2: Prediction using sampling-based visibility graph planner with A* search
    experiment_name = 'IMP_crosslinks_prediction_a_star'
    export_folder = os.path.join(cwd, 'results', experiment_name)
    threads = 4  # Number of threads for parallel processing
    print(f"\t#5.2: Predicting crosslinks with A* search (this may take a while...) → {export_folder}")

    start_time = time.time()
    prediction = structure.predict_crosslinks(
        pcd,
        from_residues,
        to_residues,
        min_distance,
        max_distance,
        'pseudoDSBU',
        direct_path=False,
        num_processes=threads
    )
    elapsed_seconds = int(time.time() - start_time)
    hours, minutes, seconds = elapsed_seconds // 3600, (elapsed_seconds % 3600) // 60, elapsed_seconds % 60
    print(f"\tA* search prediction completed in {hours}h {minutes}m {seconds}s")

    prediction.export_for_chimerax(pcd, export_folder, experiment_name)
    print("#5: Done\n")
