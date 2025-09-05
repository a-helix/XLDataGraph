# XLDataGraph

**XLDataGraph** is a feature-rich Python library for processing, filtering, comparing, and visualizing crosslinking mass spectrometry results. Built for advanced structural biology and proteomics workflows, it offers seamless integration of sequence, domain, and structural data, supporting publication-ready visualizations and network analyses.

***

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Supported File Types](#supported-file-types)
- [Quick Start Example](#quick-start-example)
- [API Overview](#api-overview)
    - [Data Loading \& Filtering](#data-loading--filtering)
    - [Visualization](#visualization)
    - [Structural Prediction](#structural-prediction)
- [Detailed Method Documentation](#detailed-method-documentation)
- [Examples](#examples)
- [Contributing](#contributing)
- [License](#license)
- [Contact](#contact)
- [Project Status](#project-status)

***

## Features

- **Crosslink data**: Imports MeroX `.zhrm` results and extracts crosslink site data.
- **Metadata**: Includes support for domain (`.dmn`), FASTA (`.fasta`), protein chain (`.pcd`), and structure (`.cif`, `.pdb`) files.
- **Flexible filtering**: By score, replica number, crosslink type (inter/intra/homeotypic).
- **Data merging/replica handling**: Directly support combining, blanking, and comparing datasets.
- **Visualization suite**:
    - Circos protein interaction plots with domain coloring.
    - Venn diagrams for cross-dataset overlap.
    - Gephi-compatible exports (AAI/PPI networks).
    - ChimeraX pseudo-bond distance constraints.
- **Structural predictions**: Fast CA-CA or advanced A* path-based 3D predictions for hypothetical crosslinks.
- **Scriptable** and extensible for custom workflows.

***

## Installation

```bash
pip install xldg
```


***

## Supported File Types

| Extension | Description |
| :-- | :-- |
| `.fasta` | Protein sequences (Uniprot, Araport, Custom) |
| `.dmn` | Protein domain annotation (CSV) |
| `.zhrm` | MeroX result (ZIP archive with CSV) |
| `.pcd` | Protein-chain assignment descriptor |
| `.cif`/`.pdb` | Structure files (mmCIF or PDB format) |


***

## Quick Start Example

```python
import os
from xldg.data import Path, MeroX, Domain, Fasta, CrossLink
from xldg.graphics import CircosConfig, Circos

cwd = './examples/files'
fasta = Fasta.load_data(os.path.join(cwd, 'example.fasta'), 'Custom')
domains = Domain.load_data(Path.list_given_type_files(cwd, 'dmn'))
crosslinks = MeroX.load_data(Path.list_given_type_files(cwd, 'zhrm'), 'DSBU')
combined = CrossLink.combine_all(crosslinks)
config = CircosConfig(fasta, domains)
circos = Circos(combined, config)
circos.save(os.path.join(cwd, 'results', 'circos_basic.svg'))
```


***

## API Overview

### Data Loading \& Filtering

#### Loading Data

- **Fasta.load_data(path, fasta_format, remove_parenthesis=False)**
    - Load one or more FASTA files. See below for argument details.
- **Domain.load_data(path)**
    - Load one or more domain annotation files.
- **MeroX.load_data(path, linker=None)**
    - Import one or more `.zhrm` crosslink result files.
- **ProteinChain.load_data(path)**
    - Load a protein-chain map (`.pcd`).
- **ProteinStructure.load_data(path)**
    - Load a 3D structure file (`.pdb` or `.cif`).


#### Filtering and Dataset Manipulation

- **CrossLink.filter_by_score(dataset, min_score=0, max_score=sys.maxsize)**
    - Keep only crosslinks in the specified score window.
- **CrossLink.filter_by_replica(dataset, min_replica=1, max_replica=sys.maxsize)**
    - Filter crosslinks by replica count.
- **CrossLink.remove_interprotein(dataset)**
    - Remove all interprotein crosslinks.
- **CrossLink.remove_intraprotein(dataset)**
    - Remove all intraprotein crosslinks.
- **CrossLink.remove_homeotypic(dataset)**
    - Remove all homeotypic crosslinks (same residue/peptide on both sides).
- **CrossLink.combine_all([datasets])**
    - Merge given datasets as a single dataset.
- **CrossLink.combine_replicas(dataset_list, n)**
    - Combine every `n` datasets into a multi-replicate dataset group.
- **CrossLink.blank_replica(dataset)**
    - Set all replica counts to 1 for plotting/overlap analyses.

***

### Visualization

- **CircosConfig**: Configures Circos protein plot visuals and filters.
- **Circos**: Generates and saves Circos plots.
- **VennConfig**: Configures Venn diagrams.
- **Venn2/Venn3**: 2- or 3-group overlap plots.
- **CrossLinkDataset.export_ppis_for_gephi(folder, filename, pcd)**
    - Exports protein-protein interaction graphs for Gephi visualization.
- **CrossLinkDataset.export_aais_for_gephi(folder, filename, pcd)**
    - Exports residue–residue networks.
- **CrossLinkDataset.export_for_chimerax(...)**
    - Exports `.pb` pseudo-bond files for ChimeraX (see below for expanded argument docs).

***

### Structural Prediction

- **ProteinStructureDataset.predict_crosslinks(...)**
    - Predicts possible crosslinks based on atomic-residue coordinates; can use direct or sampled pathfinding.

***

## Detailed Method Documentation

### Data Loading and Filtering

#### `Fasta.load_data(path, fasta_format, remove_parenthesis=False)`

- **path:** `str` or `list` of `str`. Filepaths to FASTA files.
- **fasta_format:** `str`. `"Uniprot"`, `"Araport11"`, or `"Custom"` (header parser).
- **remove_parenthesis:** `bool` (optional). Remove parentheses content from headers.
- **returns:** `FastaDataset`.


#### `Domain.load_data(path)`

- **path:** `str` or `list` of `str`. Path(s) to `.dmn` file(s).
- **returns:** `DomainDataset`.


#### `MeroX.load_data(path, linker=None)`

- **path:** `str` or `list` of `str`. Path(s) to `.zhrm` zipped result file(s).
- **linker:** `str` (optional). Linker type, for annotation purposes.
- **returns:** One `CrossLinkDataset` for single file; `list` of datasets for multiple files.


#### `CrossLink.filter_by_score(dataset, min_score=0, max_score=sys.maxsize)`

- **dataset:** `CrossLinkDataset` or `list`. Input data.
- **min_score / max_score:** `int`. Score window.
- **returns:** Filtered dataset.


#### `CrossLink.filter_by_replica(dataset, min_replica=1, max_replica=sys.maxsize)`

- **dataset:** `CrossLinkDataset` or `list`. Input data.
- **min_replica / max_replica:** `int`. Allowed replica count window.
- **returns:** Filtered dataset.


#### `CrossLink.remove_interprotein(dataset)`

- **dataset:** `CrossLinkDataset` or `list`. Dataset(s) to filter.
- **returns:** Dataset with **all crosslinks between distinct proteins removed** (only intraprotein and homeotypic remain).[^2][^1]
- **Effect:** Used to focus on internal protein organization, e.g., in Circos plots.


#### `CrossLink.remove_intraprotein(dataset)`

- **dataset:** `CrossLinkDataset` or `list`. Dataset(s) to filter.
- **returns:** Dataset with **all intra-protein crosslinks removed** (only interprotein and homeotypic links remain).[^1][^2]
- **Effect:** Good for focusing on protein–protein interactions (interactions between different proteins).


#### `CrossLink.remove_homeotypic(dataset)`

- **dataset:** `CrossLinkDataset` or `list`. Dataset(s) to filter.
- **returns:** Dataset with **all homeotypic crosslinks removed** (where both sites correspond to the same residue or peptide).[^2][^1]
- **Effect:** Streamlines network/structural analyses by removing redundancy (self/self crosslinks).


#### `CrossLink.combine_all(datasets)`

- **datasets:** `list` of `CrossLinkDataset`. All datasets to merge.
- **returns:** Combined `CrossLinkDataset`.


#### `CrossLink.combine_replicas(dataset_list, n)`

- **dataset_list:** List of datasets (e.g., by replicate).
- **n:** Number per group (e.g., `n=3` for three-replicate overlays).
- **returns:** List of merged multi-replicate datasets.


#### `CrossLink.blank_replica(dataset)`

- **dataset:** `CrossLinkDataset` or `list`. All replica counts set to 1.
- **returns:** Dataset(s) for plotting or overlap comparison.

***

### Visualization Methods

#### `CircosConfig(fasta, domains=None, ...)`

- **fasta:** FastaDataset.
- **domains:** DomainDataset (optional).
- **legend, title:** `str` (optional).
- **figsize:** `(float, float)`, e.g. `(9, 9)` for image size.
- **label_interval, space_between_sectors, font sizes, color overrides:** See docstrings/defaults in code for advanced settings.
- **returns:** CircosConfig object.


#### `Circos(crosslinks, config)`

- **crosslinks:** CrossLinkDataset.
- **config:** CircosConfig.
- **.save(path):** Renders and saves Circos plot.


#### `VennConfig(label_1, label_2, label_3=None, title=None, ...)`

- **label_1, label_2, label_3:** `str`. Categories for Venn sets (up to 3).
- **title:** `str` (optional). Plot title and additional options for color and font size.
- **returns:** VennConfig.


#### `Venn2(dataset1, dataset2, config)`

- **dataset1, dataset2:** CrossLinkDataset. Sets to compare.
- **config:** VennConfig.
- **.save(path):** Save the Venn plot.


#### `Venn3(dataset1, dataset2, dataset3, config)`

- **Add third dataset. Otherwise, usage as in Venn2.**


#### `CrossLinkDataset.export_ppis_for_gephi(folder, filename, pcd)`

- **folder:** `str`. Output folder.
- **filename:** `str`. Output `.gexf` filename.
- **pcd:** ProteinChainDataset.
- **returns:** None. Writes file.


#### `CrossLinkDataset.export_aais_for_gephi(folder, filename, pcd)`

- **Same as above, but for residue-residue/AAI level network.**


#### `CrossLinkDataset.export_for_chimerax(pcd, folder, filename, diameter=0.2, dashes=1, color_valid_distance='#48cae4', color_invalid_outsider='#d62828', protein_structure=None, min_distance=0, max_distance=sys.maxsize, atom_type='CA')`

- **pcd:** ProteinChainDataset.
- **folder:** str.
- **filename:** str.
- **diameter:** float. Bond diameter for visualization.
- **dashes:** int. Style parameter for ChimeraX.
- **color_valid_distance:** str. For valid range links.
- **color_invalid_outsider:** str. For out-of-range links.
- **protein_structure:** ProteinStructureDataset (optional). Used for distance validation.
- **min_distance/max_distance:** float. Site distance boundaries.
- **atom_type:** str. (Usually "CA").
- **returns:** None. Writes one or more `*.pb` files for ChimeraX.

***

### Structure-Based Crosslink Prediction

#### `ProteinStructureDataset.predict_crosslinks(pcd, residues_1, residues_2, min_length=1.0, max_length=sys.maxsize, linker=None, atom_type='CA', direct_path=True, radius=1.925, node_multiplier=100, num_processes=1)`

- **pcd:** ProteinChainDataset.
- **residues_1/residues_2:** str. Residue selectors, e.g. `{K` for N-term lysine, `K` for all lysines.
- **min_length:** float (angstrom). Minimum allowed CA-CA distance.
- **max_length:** float (angstrom). Maximum allowed distance.
- **linker:** str. Linker identifier.
- **atom_type:** str (default `"CA"`). On which atom dummy links will be modeled.
- **direct_path:** bool. If `True`, use direct CA-CA Euclidean distance; if `False`, use A* path sampling (models obstacles, slow).
- **radius:** float. Excluded-volume radius.
- **node_multiplier:** int. Controls sampling density if A* search used.
- **num_processes:** CPU count for parallel calculation (A* only).
- **returns:** CrossLinkDataset of predicted sites.

***

## Examples

See `examples/` directory for:

- Circos plotting
- Crosslink prediction and ChimeraX export
- Gephi (network) and Venn visualization
- Combined, advanced dataset filtering and merging

***

## Contributing

Issues and pull requests welcome! See [GitHub Issues](https://github.com/a-helix/XLDataGraph/issues).

***

## License

This project is licensed under the GNU GPLv3.

***

## Contact

- GitHub: [@a-helix](https://github.com/a-helix)

***

## Project Status

XLDataGraph is actively developed — please see the repository for the latest features, bugfixes, and documentation.
