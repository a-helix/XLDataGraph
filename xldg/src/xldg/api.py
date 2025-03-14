import os
import time
import requests
import copy
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patheffects as path_effects
import numpy as np
from typing import  List, Set, Dict
from xldg import xl


class PPI: 
# Protein-protein interaction (PPI)
    def __init__(self, prot_1: str, prot_2: str, string_id_1: str, string_id_2: str, uniprot_id_prot_1: str, string_score: str):
        self.prot_1 = prot_1
        self.prot_2 = prot_2

        self.string_id_1 = string_id_1
        self.string_id_2 = string_id_2

        self.uniprot_id_1 = uniprot_id_prot_1
        self.uniprot_id_2 = ""

        self.string_score = float(string_score)

    def __str__(self):
        return f'{self.uniprot_id_1}\t{self.uniprot_id_2}\t{self.prot_1}\t{self.prot_2}\t{self.score}\n'

    def __eq__(self, other):
        if not isinstance(other, PPI):
            return NotImplemented
        # Ensure the order of proteins doesn't matter for equality
        return (
            {self.prot_1, self.prot_2, self.uniprot_id_1, self.uniprot_id_2, self.string_id_1, self.string_id_2} == {other.prot_1, other.prot_2, other.uniprot_id_1, other.uniprot_id_2, other.string_id_1,other.string_id_2}
        )

    def __hash__(self):
        # Hash based on an unordered tuple of proteins and the score
        return hash((self.prot_1, self.prot_2, self.uniprot_id_1, self.uniprot_id_2, self.string_id_1, self.string_id_2))

    def _get_uniprot_id_by_string_id(self, string_id: str) -> str:
        # Make the request
        params = {
            "query": string_id,
            "fields": [
            "accession",
            "protein_name",
            "cc_function",
            "ft_binding"
            ],
            "sort": "accession desc",
            "size": "50"
        }
        headers = {
          "accept": "application/json"
        }
        base_url = "https://rest.uniprot.org/uniprotkb/search"

        response = requests.get(base_url, headers=headers, params=params)
        if response.status_code == 200:
            data = response.json()
            if "results" in data:
                for result in data["results"]:
                    primary_accession = result.get("primaryAccession")
                    return primary_accession
        else:
            print(f"Error: {response.status_code} - {response.text}")
    
    def update_uniprot_id_by_string_id(self):
        self.uniprot_id_2 = self._get_uniprot_id_by_string_id(self.string_id_2)
        
        if self.uniprot_id_2 == None:
            self.uniprot_id_2 = self._get_uniprot_id_by_string_id(f'{self.prot_2} Homo sapiens (Human)') # TODO
        if self.uniprot_id_2 == None:
            self.uniprot_id_2 = self._get_uniprot_id_by_string_id(f'{self.prot_2}_HUMAN')
        if self.uniprot_id_2 == None:
            stripped_id = self.string_id_2.replace("9606.", "")
            self.uniprot_id_2 = self._get_uniprot_id_by_string_id(stripped_id)
        if self.uniprot_id_2 == None:
            print(f'ERROR! Unable to find Uniprot ID for STRING ID {self.string_id_2}.')


def save_raw_string_PPIs(raw_save_path: str, proteins: Set['Uniprot_Fasta'], species_ncbi_id: str, sleep: float = 0.5) -> None:
    string_api_url = "https://version-12-0.string-db.org/api"
    output_format = "tsv-no-header"
    method = "PPI_partners"

    request_url = "/".join([string_api_url, output_format, method])

    protein_uniprot_ids = [prot.uniprot_id for prot in proteins]

    response = []
    counter = 1
    for uniprot_id in protein_uniprot_ids:
        params = {
            "identifiers" : uniprot_id,  # your protein
            "species" :species_ncbi_id,  # species NCBI identifier 
            "limit" : 100000000,
            "caller_identity" : "xldg"
        }
        time.sleep(sleep)
        
        print(f"   {counter} of {len(protein_uniprot_ids)}. Running {uniprot_id} STRING database request.")
        counter += 1

        string_db_output = requests.post(request_url, data=params).text
        if "Error" in string_db_output:
            for entity in proteins:
                if uniprot_id == entity.uniprot_id:
                    print(f"ERROR!. Rejected STRING database request. {uniprot_id} from {next(iter(entity.pathway))}")
                    break
            continue

        data = string_db_output.strip().split("\n")
        buffer = []

        for row in data:
            buffer.append(f"{row}\t{uniprot_id}")
            
        response += buffer

    # Write PPIs to the specified file
    with open(raw_save_path, 'a') as file:
        for line in response:
            file.write(line + '\n')
    print(f'Raw file saved to {raw_save_path}')


def filter_raw_data_by_score(raw_data_path: str, save_path: str, filter_score: float = 0.7) -> None: 
    # TODO: Add min and max score filter
    raw_data = []
    with open(raw_data_path, 'r') as file:
        for line in file:
            raw_data.append(line)

    all_interactors = []
    print("Identifying all PPIs.")

    for line in raw_data:
        buffer = []
        l = line.strip().split("\t")
        string_id_1 = l[0]
        string_id_2 = l[1]
        query_name = l[2]
        partner_name = l[3]
        combined_score = l[5]
        uniprot_id_prot_1 = l[-1]

        buffer.append(query_name)
        buffer.append(partner_name)
        buffer = sorted(buffer)

        all_interactors.append(PPI(buffer[0], buffer[1], string_id_1, string_id_2, uniprot_id_prot_1, combined_score))

    buffer = set()

    print(f"Filtering PPIs using score {filter_score}.")
    for PPI in all_interactors:
        if PPI.score >= filter_score:
            buffer.add(PPI)
    
    print("Updating all Uiprot IDs. This might take some time.")
    len_buffer = len(buffer)
    counter = 0

    for unit in buffer:
        unit.update_uniprot_id_by_string_id()

        counter += 1
        if counter % 100 == 0:
            print(f"\t{counter} of {len_buffer} PPIs were updated.")

    with open(save_path, 'w') as file:
        for unit in buffer:
            file.write(str(unit))
    print(f'Processed data saved to {save_path}')

class Edge:
    def __init__(self, uniprot_id_1: str, uniprot_id_2: str, protein_1: str, protein_2: str): 
        self.uniprot_id_1 = uniprot_id_1
        self.uniprot_id_2 = uniprot_id_2

        self.protein_1 = protein_1
        self.protein_2 = protein_2

        self.count = 0
        self.protein_1_necro_path = set()
        self.protein_2_necro_path = set()
        self.common_necro_pathways = set()

    def __hash__(self):
        return hash(tuple(sorted((self.protein_1, self.protein_2, self.uniprot_id_1, self.uniprot_id_2))))

    def __eq__(self, other):
        return (sorted((self.protein_1, self.protein_2, self.uniprot_id_1, self.uniprot_id_2)) == 
                sorted((other.protein_1, other.protein_2, other.uniprot_id_1, other.uniprot_id_2)))

    def __str__(self):
        return f'{self.uniprot_id_1}\t{self.uniprot_id_2}\t{self.protein_1}\t{self.protein_2}\t{self.count}\t{self.common_necro_pathways}'

    @classmethod
    def from_tsv(cls, line: str) -> 'Edge':
        # Split TSV line into fields
        fields = line.split("\t")
        
        # Extract mandatory fields
        uniprot_id_1 = fields[0]
        uniprot_id_2 = fields[1]
        protein_1 = fields[2]
        protein_2 = fields[3]
        
        # Create an instance
        edge = cls(uniprot_id_1, uniprot_id_2, protein_1, protein_2)
        
        # Parse count as integer
        edge.count = int(fields[4])

        # Parse necroptosis pathways as sets
        necro_pathways = eval(fields[5])  # Convert string representation of a set to an actual set
        edge.protein_1_necro_path = necro_pathways
        edge.protein_2_necro_path = necro_pathways
        edge.common_necro_pathways = edge.protein_1_necro_path.intersection(edge.protein_2_necro_path)
        
        return edge

def generate_edges(
    path_to_edges: str, 
    input_fasta_set: Set['Uniprot_Fasta'], 
    save_path: str = None
) -> List['Edge']:
    fasta_set = copy.deepcopy(input_fasta_set)

    sorted_fasta_set = sorted(input_fasta_set, key=lambda x: (
        len(x.pathway), 
        x.gene, 
        x.uniprot_id  # Third, stable sorting key
    ))
    
    # Create dictionary with sorted, stable keys
    fasta_dict = {fasta.uniprot_id: sorted(fasta.pathway) for fasta in sorted_fasta_set}

    print(f'Generating edges from {path_to_edges}')
    edges = []

    # Load edges from the file
    try:
        with open(path_to_edges, "r") as file:
            for line in file:
                data = line.strip().split('\t')
                if len(data) < 4:  # Ensure there's enough data to create an edge
                    print(f"Skipping malformed line: {line}")
                    continue
                
                prot_name_1 = ''
                uniprot_id_1 = ''

                prot_name_2 = ''
                uniprot_id_2 = ''

                for fasta in input_fasta_set:
                    if fasta.uniprot_id == data[0] and fasta.gene == data[2]:
                        prot_name_1 = fasta.gene
                        uniprot_id_1 = fasta.uniprot_id

                    elif fasta.uniprot_id == data[1] and fasta.gene == data[3]:
                        prot_name_2 = fasta.gene
                        uniprot_id_2 = fasta.uniprot_id

                    if prot_name_1 != '' and prot_name_2 != '':
                        break

                    if uniprot_id_2 == 'None':
                        break

                edge = Edge(uniprot_id_1, uniprot_id_2, prot_name_1, prot_name_2)
                edges.append(edge)

    except FileNotFoundError:
        print(f"Error: File {path_to_edges} not found.")
        return []

    except Exception as e:
        print(f"Error reading file {path_to_edges}: {e}")
        return []

    updated_edges = []

    for edge in edges:
        uniprot_1 = edge.uniprot_id_1
        uniprot_2 = edge.uniprot_id_2

        if uniprot_1 not in fasta_dict or uniprot_2 not in fasta_dict:
            continue

        path_1 = fasta_dict.get(uniprot_1, None)
        path_2 = fasta_dict.get(uniprot_2, None)

        if path_1 is None or path_2 is None:
            continue

        edge.protein_1_necro_path.update(path_1)
        edge.protein_2_necro_path.update(path_2)

        # Compute common pathways
        edge.common_necro_pathways = set(sorted(set(path_1) & set(path_2)))
        edge.count = len(edge.common_necro_pathways)

        # Only add edges with meaningful pathways
        if edge.count > 0:
            updated_edges.append(edge)
    
    # Save edges to file (if save_path is provided)
    updated_edges = sorted(updated_edges, key=lambda edge: edge.count, reverse=True)
    if save_path:
        print(f'Saving edges to {save_path}')
        try:
            with open(save_path, 'w') as file:
                for edge in updated_edges:
                    file.write(f'{edge.uniprot_id_1}\t{edge.uniprot_id_2}\t{edge.protein_1}\t{edge.protein_2}\t{edge.count}\t{str(edge.common_necro_pathways)}\n')               

        except Exception as e:
            print(f"Error saving edges to {save_path}: {e}")
    return updated_edges

def download_edges_from_tsv_file(path: str) -> List['Edge']:
    edges = []
    
    with open(path, "r") as file:
        for line in file:
            # Skip empty lines or headers if present
            if line.strip():
                edges.append(Edge.from_tsv(line.strip()))
    
    return edges

def download_alphafold_pdb(uniprot_id: str, prot_name: str, save_dir: str = "AlphaFold3_structures") -> str:
    base_url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
    
    # Create save directory if it doesn't exist
    os.makedirs(save_dir, exist_ok=True)
    
    # Path to save the PDB file
    pdb_file_path = os.path.join(save_dir, f"{prot_name}({uniprot_id}).pdb")
    
    try:
        # Download the PDB file
        response = requests.get(base_url)
        response.raise_for_status()  # Raise an exception for HTTP errors
        
        # Save the file
        with open(pdb_file_path, "wb") as pdb_file:
            pdb_file.write(response.content)
        
        print(f"Downloaded: {pdb_file_path}")
        return pdb_file_path
    except requests.exceptions.HTTPError as http_err:
        print(f"\t{http_err}")
        return ""
    except Exception as err:
        print(f"Unknown error: {err}")
        return ""

def download_protein_pdbs_of_edges(edges: List['Edge'], input_fasta_set: Set['Uniprot_Fasta'], min_path: int = 2, save_dir: str = "AlphaFold3_structures") -> List[str]:
    saved_files = []

    for edge in edges:
        if len(edge.common_necro_pathways) < min_path:
            continue

        prot_1_file_path = os.path.join(save_dir, f"{edge.protein_1}({edge.uniprot_id_1}).pdb")
        if os.path.exists(prot_1_file_path) is False:
            path = download_alphafold_pdb(edge.uniprot_id_1, edge.protein_1, save_dir)
            if path != "":
                saved_files.append(path)
        else:
            saved_files.append(prot_1_file_path)
            

        prot_2_file_path = os.path.join(save_dir, f"{edge.protein_2}({edge.uniprot_id_2}).pdb")

        if os.path.exists(prot_2_file_path) is False:
            path = download_alphafold_pdb(edge.uniprot_id_2, edge.protein_2, save_dir)
            if path != "":
                saved_files.append(path)
        else:
            saved_files.append(prot_2_file_path)
            
    print(f"PDBs were saved to \"{save_dir}\"")
    return saved_files

class pLDDT_sequence:
    def __init__(self, start: int = 0, end: int = 0, pLDDT_category: str = ''):
        self.start = start
        self.end = end
        self.pLDDT_category = pLDDT_category

    def __str__(self):
        string = ''
        return f'{self.start}\t{self.end}\t{self.pLDDT_category}'

    def __eq__(self, other):
        if not isinstance(other, pLDDT_sequence):
            return NotImplemented
        return (
            self.start == other.start and
            self.end == other.end and
            self.pLDDT_category == other.pLDDT_category
        )

    def __hash__(self):
        """Generate a hash for the pLDDT_sequence."""
        return hash((self.start, self.end, self.pLDDT_category))

class Protein_Disorder_Map:
    def __init__(self, pdb_path: str, scale_bar = False):
        self.path = pdb_path
        self.name = ''
        self.uniprot_id = ''
        self.disordered_aa = 0
        self.ordered_aa = 0
        self.disordered_percentage = 0
        self.pLDDT_sequences = []

        if scale_bar == False:
            self._extract_prot_name_and_uniprot_id()
            self._map_disorder_from_pdb()
            self.disordered_percentage = self.disordered_aa / (self.disordered_aa + self.ordered_aa)
        else:
             self._initialise_scale_bar()
            
    def __str__(self):
        joined_pLDDT_sequences = "\n".join(str(obj) for obj in self.pLDDT_sequences)
        return f'{self.path}\n{self.name}\n{self.uniprot_id}\n{self.disordered_aa}\n{self.ordered_aa}\n{self.disordered_percentage}\n{joined_pLDDT_sequences}\n'

    def __eq__(self, other):
        if not isinstance(other, Protein_Disorder_Map):
            return NotImplemented
        return (
            self.path == other.path and
            self.name == other.name and
            self.uniprot_id == other.uniprot_id and
            self.disordered_aa == other.disordered_aa and
            self.ordered_aa == other.ordered_aa and
            self.pLDDT_sequences == other.pLDDT_sequences and
            abs(self.disordered_percentage - other.disordered_percentage) < 1e-6
        )

    def __hash__(self):
        """Generate a hash for the ProteinDisorderMap."""
        return hash((
            self.path,
            self.name,
            self.uniprot_id,
            self.disordered_aa,
            self.ordered_aa,
            tuple(self.pLDDT_sequences),  # Convert the list of sequences to a tuple for hashing
            round(self.disordered_percentage, 6)  # Avoid floating-point issues in hash
        ))

    def _extract_prot_name_and_uniprot_id(self) -> None:
        file_name_without_ext = os.path.splitext(os.path.basename(self.path))[0]
        splited_name = file_name_without_ext.replace(')', '').split('(')
        self.name = splited_name[0]
        self.uniprot_id = splited_name[1]

    def _identify_pLDDT_category_by_score(self, score: float) -> str:
        pLDDT_categorie = {
            "VERY_LOW": 50.0,
            "LOW": 70.0,
            "HIGH": 90.0,
            "VERY_HIGH": 100.0
        }

        for category, limit in pLDDT_categorie.items():
            if score < limit:
                return category

    def _map_disorder_from_pdb(self):
        with open(self.path, "r") as file:
            buffer = []
            for line in file:
                splited_line = line.strip().split()
                if splited_line[0] == "ATOM" and splited_line[2] == "CA":
                    buffer.append(splited_line)

            counter = 0 
            sequence = pLDDT_sequence()
            pLDTT_score = 0
            category = ''
            residue_number = 0

            for splited_line in buffer:
                counter += 1
                pLDTT_score = float(splited_line[-2])
                category = self._identify_pLDDT_category_by_score(pLDTT_score)

                if "." in splited_line[5]:
                    residue_number = int(splited_line[4].replace("A", ""))
                else:
                    residue_number = int(splited_line[5])

                if pLDTT_score < 70.0:
                    self.disordered_aa += 1
                else:
                    self.ordered_aa += 1
                   
                if sequence.pLDDT_category == '':
                    sequence.pLDDT_category = category
                    sequence.start = residue_number
                    continue

                if sequence.pLDDT_category != category:
                    sequence.end = residue_number
                    self.pLDDT_sequences.append(sequence)
                    sequence = pLDDT_sequence()
                    category = ''
                    continue

                if counter == len(buffer):
                    sequence.end = residue_number
                    self.pLDDT_sequences.append(sequence)

    def _initialise_scale_bar(self) -> None:
        self.name = "Scale Bar"
        self.uniprot_id = "250 AA"
        self.pLDDT_sequences.append(pLDDT_sequence(1, 250, "SCALE_BAR"))


def extract_disorder_information_from_pdb_files(path_to_files: List[str]) -> List['Protein_Disorder_Map']:
    buffer = []
    for path in path_to_files:
        buffer.append(Protein_Disorder_Map(path))
    return buffer

def create_protein_info_table(protein_maps: List['Protein_Disorder_Map'], output_file: str = "protein_disorder_table.svg") -> None:
    # Sort input by disorder percentage (descending) and remove duplicates
    protein_maps = sorted(set(protein_maps), key=lambda x: x.disordered_percentage, reverse=True)

    # Prepare table data
    cell_text = []
    for idx, protein in enumerate(protein_maps, start=1):
        cell_text.append([
            idx,
            protein.name,
            protein.uniprot_id,
            f"{protein.disordered_percentage * 100:.2f}",
            protein.pLDDT_sequences[-1].end
        ])

    # Calculate the number of rows in the table (including header)
    n_rows = len(protein_maps) + 1  # Including header row

    # Adjust the figure size based on the number of rows in the table
    fig_height = max(3, n_rows * 0.4)  # Minimum height of 3 inches for small tables
    fig_width = 16  # Adjusted width for the figure

    # Create figure
    fig = plt.figure(figsize=(fig_width, fig_height))

    # Create a gridspec with padding
    gs = fig.add_gridspec(1, 2, width_ratios=[1.2, 1], wspace=0.1)
    ax = fig.add_subplot(gs[0])  # Table

    # Remove axes
    ax.set_axis_off()
    # Create table with adjusted position
    table = ax.table(
        cellText=cell_text,
        colLabels=['#', 'Gene', 'UniProt ID', 'Disorder (%)', 'Size (AA)'],
        loc='center',
        cellLoc='center',
        bbox=[0, 0.00, 1.0, 1.0],
        colWidths=[0.06, 0.27, 0.25, 0.25, 0.17]  # Column widths adjusted
    )

    # Style the table
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    
    # Calculate row heights for alignment
    row_height = 1 / n_rows
    table.scale(1, row_height * 2.5)

    # Style cells
    for i in range(n_rows):
        for j in range(5):
            cell = table[(i, j)]
            if i == 0:  # Header
                cell.set_facecolor('#e0e0e0')
                cell.set_text_props(weight='bold')
            else:  # Data rows
                cell.set_facecolor('#ffffff')
            cell.set_edgecolor('black')

    # Save the figure with high resolution
    plt.savefig(output_file, dpi=600, bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    print('START')
    fasta_folder_path = r"path_to_fasta_folder"
    all_fasta_files_locations = list_fasta_type_files_in_folder(fasta_folder_path)
    all_fasta_files_content = read_all_uniprot_fasta_from_path_list(all_fasta_files_locations)

    necroptome_pathways = get_all_necroptome_pathways(all_fasta_files_content)

    raw_tsv_path = os.path.join(fasta_folder_path, r'raw_PPIs.tsv')
    # save_raw_string_PPIs(raw_tsv_path, all_fasta_files_content)
    
    processed_tsv_path_07 = os.path.join(fasta_folder_path, r'filtered_PPIs_07.tsv')
    # filter_raw_data_by_score(raw_tsv_path, processed_tsv_path_07, 0.7)

    processed_tsv_path_09 = os.path.join(fasta_folder_path, r'filtered_PPIs_09.tsv')
    # filter_raw_data_by_score(raw_tsv_path, processed_tsv_path_09, 0.9)

    # Colors for heatmap
    colors = [
        "#FF5733", "#33FF57", "#3357FF", "#FF33A6", "#F3FF33",
        "#FF8C33", "#33FFF6", "#F633FF", "#33FF8C", "#FF3333",
        "#33D4FF", "#8C33FF", "#FF8333", "#33FF44"
    ]

    necro_path = list(necroptome_pathways)

    # 0.7 STRING score
    save_edges_07 = os.path.join(fasta_folder_path, r'filtered_PPIs_with_RCD_pathways_07.tsv')
    # edges07 = generate_edges(processed_tsv_path_07, all_fasta_files_content, save_edges_07)
    edges07 = download_edges_from_tsv_file(save_edges_07)
    edges07 = [edge for edge in edges07 if edge.count > 1]

    # 0.9 STRING score
    save_edges_09 = os.path.join(fasta_folder_path, r'filtered_PPIs_with_RCD_pathways_09.tsv')
    # edges09 = generate_edges(processed_tsv_path_09, all_fasta_files_content, save_edges_09)
    edges09 = download_edges_from_tsv_file(save_edges_09)
    edges09 = [edge for edge in edges09 if edge.count > 1]

    ### Infographics for PPIs with >1 pathway
    folder = os.path.join(os.getcwd(), 'Two_and_more_RCD_pathways_PPIs')
    # 0.7 STRING score
    edges07_pdb_path_list = download_protein_pdbs_of_edges(edges07, all_fasta_files_content, min_path=2)
    edges07_disorder_info = extract_disorder_information_from_pdb_files(edges07_pdb_path_list)
    create_protein_info_table(edges07_disorder_info, f"{folder}\\Figure_S3_proteins_disorder_table_07.svg")
    create_proteins_alphafold_disorder_scheme(edges07_disorder_info, f"{folder}\\Figure_S1_alphafold3_disorder_summary_07.svg")

    plot_PPIs_heatmap_in_sections(colors, necro_path, edges07, rows_in_section=100, save_path=f"{folder}\\Figure_S5_PPIs_map_07.svg")
    plot_PPI_frequency_heatmap(edges07, necro_path, save_path=f"{folder}\\Figure_S7_RCD_pathways_summary_07.svg")

    # 0.9 STRING score
    edges09_pdb_path_list = download_protein_pdbs_of_edges(edges09, all_fasta_files_content, min_path=2)
    edges09_disorder_info = extract_disorder_information_from_pdb_files(edges09_pdb_path_list)
    create_protein_info_table(edges09_disorder_info, f"{folder}\\Figure_S4_proteins_disorder_table_09.svg")
    create_proteins_alphafold_disorder_scheme(edges09_disorder_info, f"{folder}\\Figure_S2_alphafold3_disorder_summary_09.svg")

    plot_PPIs_heatmap_in_sections(colors, necro_path, edges09, rows_in_section=100,  save_path=f"{folder}\\Figure_S6_PPIs_map_09.svg")
    plot_PPI_frequency_heatmap(edges09, necro_path, save_path=f"{folder}\\Figure_S8_RCD_pathways_summary_09.svg")

    ### Infographics for PPIs with >2 pathways
    folder = os.path.join(os.getcwd(), 'Three_and_more_RCD_pathways_PPIs')
    # 0.7 STRING score
    edges07_pdb_path_list = download_protein_pdbs_of_edges(edges07, all_fasta_files_content, min_path=3)
    edges07_disorder_info = extract_disorder_information_from_pdb_files(edges07_pdb_path_list)
    create_protein_info_table(edges07_disorder_info, f"{folder}\\Figure_S11_proteins_disorder_table_07.svg")
    create_proteins_alphafold_disorder_scheme(edges07_disorder_info, f"{folder}\\Figure_S9_alphafold3_disorder_summary_07.svg")

    edges07 = [edge for edge in edges07 if edge.count > 2]
    plot_PPIs_heatmap_in_sections(colors, necro_path, edges07, rows_in_section=0, save_path=f"{folder}\\Figure_S13_PPIs_map_07.svg")
    plot_PPI_frequency_heatmap(edges07, necro_path, save_path=f"{folder}\\Figure_S15_RCD_pathways_summary_07.svg")

    # 0.9 STRING score
    edges09_pdb_path_list = download_protein_pdbs_of_edges(edges09, all_fasta_files_content, min_path=3)
    edges09_disorder_info = extract_disorder_information_from_pdb_files(edges09_pdb_path_list)
    create_protein_info_table(edges09_disorder_info, f"{folder}\\Figure_S12_proteins_disorder_table_09.svg")
    create_proteins_alphafold_disorder_scheme(edges09_disorder_info, f"{folder}\\Figure_S10_alphafold3_disorder_summary_09.svg")

    edges09 = [edge for edge in edges09 if edge.count > 2]
    plot_PPIs_heatmap_in_sections(colors, necro_path, edges09, rows_in_section=0,  save_path=f"{folder}\\Figure_S14_PPIs_map_09.svg")
    plot_PPI_frequency_heatmap(edges09, necro_path, save_path=f"{folder}\\Figure_S16_RCD_pathways_summary_09.svg")
    print('FINISHED')

