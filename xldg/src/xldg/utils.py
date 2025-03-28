from dataclasses import dataclass
import zipfile
import os
import io
import re
import sys
import copy
from typing import List, Tuple, Set, Union, Optional
import xml.etree.ElementTree as ET
from xml.dom import minidom

from xldg.xl import CrossLink, CrossLinkDataset


class Path:
    @staticmethod
    def list_specified_type_files_from_folder(folder_path: str, file_format: str) -> List[str]:
        """
        Retrieve a list of files with a specific file extension from a given folder.

        This method scans the specified folder and collects full paths of files 
        that end with the given file format (extension).

        Args:
            folder_path (str): The directory path to search for files.
            file_format (str): The file extension to filter (e.g., '.txt', 'mgf').

        Returns:
            List[str]: A list of full file paths matching the specified format.

        Raises:
            FileNotFoundError: If the specified folder does not exist.
        """

        Path.validate_folder_existence(folder_path)

        files = []
        for file in os.listdir(folder_path):
            if file.endswith(file_format):
                files.append(os.path.join(folder_path, file))

        return files

    @staticmethod
    def sort_filenames_by_first_integer(strings: List[str], ignore: Optional[str] = None) -> List[str]:
        """
        Sort a list of file paths based on the first integer found in their filenames.

        This method extracts the first integer from each filename and uses it as 
        the sorting key. Files without integers are placed at the end of the list.

        Args:
            strings (List[str]): List of file paths to be sorted.
            ignore (str, optional): A substring to remove from filenames before 
                                    extracting the integer. Defaults to None.

        Returns:
            List[str]: Sorted list of file paths based on the first integer found.

        Note:
            - If no integer is found in a filename, it is assigned a maximum value, 
              effectively placing it at the end of the list.
            - The sorting is based on the numerical value of the first integer.
        """

        def _extract_integer(s: str) -> int:
            file_name = os.path.basename(s)
            if ignore is not None:
                file_name = file_name.replace(ignore, '')

            match = re.search(r'(\d+)', file_name)
            return int(match.group(1)) if match else float('inf')

        return sorted(strings, key=_extract_integer)

    @staticmethod
    def validate_file_existence(path: str) -> None:
        """
        Validate that a given path exists and is a file.

        Args:
            path (str): The file path to validate.

        Raises:
            FileNotFoundError: If the path does not exist or is not a file.
        """

        if not os.path.isfile(path):
                raise FileNotFoundError(f'File not found: {path}')

    @staticmethod
    def validate_folder_existence(path: str) -> None:
        """
        Validate that a given path exists and is a directory.

        Args:
            path (str): The folder path to validate.

        Raises:
            FileNotFoundError: If the path does not exist or is not a directory.
        """

        if not os.path.isdir(path):
                raise FileNotFoundError(f'Folder not found: {path}')

    @staticmethod
    def validate_file_format(path: str, format: str) -> None:
        """
        Validate that a file has the expected file extension.

        This method checks if the file extension matches the expected format, 
        performing a case-insensitive comparison.

        Args:
            path (str): The full path of the file to validate.
            format (str): The expected file extension.

        Raises:
            ValueError: If the file extension does not match the expected format.

        Example:
            validate_file_format('document.PDF', 'pdf')  # Passes
            validate_file_format('image.jpg', 'png')     # Raises ValueError
        """

        _, extension = os.path.splitext(path)
        extension = extension.lower().lstrip('.')
        expected_format = format.lower()
        
        if extension != expected_format:
            raise ValueError(f'Invalid file format. Expected {expected_format.upper()}, '
                             f'got {extension.upper()}')

class Data:
    @staticmethod
    def read_merox_file(path: str, linker: Optional[str] = None) -> 'CrossLinkDataset':
        """
        Read and parse a MeroX result file (.zhrm) into a CrossLinkDataset.

        This method extracts cross-linking information from a compressed MeroX result file.
        It opens the 'Result.csv' file within the zip archive and creates CrossLink objects
        for each row of data.

        Args:
            path (str): Full path to the .zhrm file (MeroX result file).
            linker (str, optional): Specific linker used in the cross-linking experiment. 
                                    Defaults to None.

        Returns:
            CrossLinkDataset: A dataset containing all cross-links found in the file.

        Raises:
            ValueError: If the file format is not '.zhrm'.
            zipfile.BadZipFile: If the zip file is corrupted.
            IOError: If there are issues reading the file.

        Note:
            - Uses 'Result.csv' within the zip file as the source of cross-link information.
            - Sets the cross-link counter to 1 for all entries in the dataset.
            - Assumes a specific column order in the CSV file for parsing cross-link data.
        """

        Path.validate_file_format(path, '.zhrm')
        xls = []
        software = 'MeroX'

        with zipfile.ZipFile(path, 'r') as zip_ref:
            with zip_ref.open('Result.csv') as csv_file:
                for line in io.TextIOWrapper(csv_file, encoding='utf-8'):
                    row = line.strip().split(';')
                    xl = CrossLink(row[7], row[6], row[8], row[9], row[20], 
                            row[11], row[10], row[12], row[13], row[21],
                            row[0], software, linker)
                    xls.append(xl)

        dataset = CrossLinkDataset(xls)
        dataset.set_xls_counter_to_one()
        return dataset

    @staticmethod
    def read_all_merox_files(path_list: List[str], linker: Optional[str] = None) -> List['CrossLinkDataset']:
        """
        Read multiple MeroX result files and create a list of CrossLinkDatasets.

        This method processes a list of .zhrm files, extracting cross-link information 
        from each file and creating a separate CrossLinkDataset for each.

        Args:
            path_list (List[str]): List of full paths to .zhrm files to be processed.
            linker (str, optional): Specific linker used in the cross-linking experiments. 
                                    Defaults to None and will be passed to each file reading.

        Returns:
            List[CrossLinkDataset]: A list of CrossLinkDataset objects, 
                                    one for each processed input file.

        Raises:
            FileNotFoundError: If any of the specified files do not exist.

        Note:
            - Prints the path of each file being extracted for user feedback.
            - Uses read_merox_file() method to process individual files.
            - Validates file existence before attempting to read.
        """
        file_content = []
        for path in path_list:
            Path.validate_file_existence(path)
            print(f'Extracting: {path}')
            file_content.append(Data.read_merox_file(path, linker))   
        return file_content

    @staticmethod
    def filter_by_score(
        dataset: Union['CrossLinkDataset', List['CrossLinkDataset']], 
        min_score: int = 0, 
        max_score: int = sys.maxsize
    ) -> Union['CrossLinkDataset', List['CrossLinkDataset']]:
        """
        Filter dataset(s) by score range.
    
        Args:
            dataset (CrossLinkDataset or List[CrossLinkDataset]): 
                Single dataset or list of datasets to filter
            min_score (int, optional): Minimum score threshold. Defaults to 0.
            max_score (int, optional): Maximum score threshold. Defaults to sys.maxsize.
    
        Returns:
            Filtered dataset(s) of the same type as input
    
        Raises:
            ValueError: If max_score is less than min_score
            TypeError: If input is not a dataset or list of datasets
        """
        if max_score < min_score:
            raise ValueError('ERROR! max_score is smaller than min_score')
    
        if isinstance(dataset, CrossLinkDataset):
            dataset.filter_by_score(min_score, max_score)
            return dataset
    
        if isinstance(dataset, list):
            for data in dataset:
                if not isinstance(data, CrossLinkDataset):
                    raise TypeError(f'Expected CrossLinkDataset, got {type(data)}')
                data.filter_by_score(min_score, max_score)
            return dataset

        raise TypeError(f'Expected CrossLinkDataset or List[CrossLinkDataset], got {type(dataset)}')

    @staticmethod
    def combine_replicas(dataset: List['CrossLinkDataset'], n = 3) -> List['CrossLinkDataset']:
        """
        Combines replicas of datasets into groups of size n.

        This method takes a list of CrossLinkDataset instances and groups them into 
        combined datasets, with each group containing exactly n datasets. If the total 
        number of datasets is not divisible by n, a ValueError is raised.

        Args:
            dataset (List[CrossLinkDataset]): List of datasets to be combined
            n (int, optional): Number of datasets to combine in each group. Defaults to 3.

        Returns:
            List[CrossLinkDataset]: A list of combined datasets, where each combined 
            dataset is created from n original datasets.

        Raises:
            ValueError: If the total number of datasets is not evenly divisible by n.
        """

        combined_dataset = []
        buffer = []
    
        if ((len(dataset) % n) != 0):
            raise ValueError(f'Dataset size {len(dataset)} is not mutiple to n={n}')
    
        for data in dataset:
            if (len(buffer) == n):
                combined_dataset.append(CrossLinkDataset.combine_datasets(buffer))
                buffer.clear()
            buffer.append(data)
        
        combined_dataset.append(CrossLinkDataset.combine_datasets(buffer))

        return combined_dataset

    @staticmethod
    def combine_all_datasets(dataset_list: List['CrossLinkDataset']) -> 'CrossLinkDataset':
        """
        Combines multiple CrossLinkDataset instances into a single dataset.

        This method iteratively combines all datasets in the input list by 
        using the '+=' operator of the CrossLinkDataset class.

        Args:
            dataset_list (List[CrossLinkDataset]): List of datasets to be combined

        Returns:
            CrossLinkDataset: A single dataset containing all data from the input list
        """
        combined_xls = None
        for dataset in dataset_list:
            if combined_xls is None:
                combined_xls = dataset
            else:
                combined_xls += dataset
        return combined_xls

    @staticmethod
    def generate_custom_list_with_int_ranges(*diapason: Tuple[int, int]) -> List[int]:
        """
        Generates a list of integers by including full ranges from multiple start-end pairs.

        This method creates a list of integers by expanding multiple integer ranges. 
        Each range is specified as a (start, end) tuple, and the method includes 
        all integers from start to end (inclusive).

        Args:
            *diapason (Tuple[int, int]): Variable number of (start, end) tuples 
            representing integer ranges

        Returns:
            List[int]: A flattened list of integers covering all specified ranges

        Raises:
            ValueError: If any start value is greater than its corresponding end value
        """

        custom_list = []

        for start, end in diapason:
            if start > end:
                raise ValueError("ERROR! start value is greater than end value")
            custom_list.extend(range(start, end + 1))

        return custom_list

    @staticmethod
    def combine_selected_datasets(dataset_list: List['CrossLinkDataset'], indexes: List[int]) -> 'CrossLinkDataset':
        """
        Combines datasets from the input list using specified indexes.

        This method allows selecting and combining specific datasets from a list 
        by their index positions.

        Args:
            dataset_list (List[CrossLinkDataset]): List of datasets to select from
            indexes (List[int]): List of indexes indicating which datasets to combine

        Returns:
            CrossLinkDataset: A single dataset created by combining the selected datasets

        Raises:
            IndexError: If any provided index is out of range of the dataset list
        """

        biggest_index = max(indexes)
        if biggest_index >= len(dataset_list):
            raise IndexError(f"ERROR! index {biggest_index} out of given dataset_list range")

        buffer = []
        for x in indexes:
            buffer.append(dataset_list[x])

        return CrossLinkDataset.combine_datasets(buffer) 

class ProteinChainDataset:
    def __init__(self, pcd_file_path: str):
        self.path = pcd_file_path
        self.pcds = {}
        self._assign_pcds(pcd_file_path)

    def _assign_pcds(self, path_to_pcd_file: str) -> None:
        Path.validate_file_existence(path_to_pcd_file)
        Path.validate_file_format(path_to_pcd_file, '.pcd')

        with open(path_to_pcd_file, 'r') as file:
            for line in file:
                splited_line = line.replace('\n', '').split(',')

                if len(splited_line) < 2:
                    raise ValueError(f'Wrong data format in {path_to_pcd_file}')

                self.pcds[splited_line[0]] = splited_line[1:]


    def __len__(self):
        return len(self.pcds)

    def __iter__(self):
        return iter(self.pcds.items())

    def __getitem__(self, key):
        return self.pcds[key]

    def __next__(self):
        return next(iter(self.pcds.items()))

@dataclass
class Atom:
    number: int
    residue: str  # Residue identifier (number + insertion code if present)
    type: str     # (e.g., 'N', 'CA')
    chain: str    
    x: float
    y: float
    z: float

class ProteinStructure:
    def __init__(self, file_path: str):
        self.file_path = file_path
        self.atoms: List[Atom] = []
        content = self._read_file()
        if file_path.lower().endswith('.pdb'):
            self._parse_pdb(content)
        elif file_path.lower().endswith('.cif'):
            self._parse_cif(content)
        else:
            raise ValueError("Unsupported file format. Only .pdb and .cif are supported.")

    def _read_file(self) -> str:
        with open(self.file_path, 'r') as f:
            return f.read()

    def _parse_pdb(self, content: str):
        for line in content.split('\n'):
            if line.startswith(('ATOM  ', 'HETATM')):
                try:
                    # Extract fields using fixed column positions
                    atom_number = int(line[6:11].strip())
                    atom_type = line[12:16].strip()
                    residue_name = line[17:20].strip()
                    chain = line[21].strip()
                    residue_number = line[22:26].strip()  # Includes insertion code
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    self.atoms.append(Atom(
                        number=atom_number,
                        residue=residue_number,
                        type=atom_type,
                        chain=chain,
                        x=x,
                        y=y,
                        z=z
                    ))
                except (ValueError, IndexError):
                    continue  # Skip invalid lines

    def _parse_cif(self, content: str):
        lines = content.split('\n')
        in_atom_site_loop = False
        headers = []
        id_idx = label_atom_id_idx = label_comp_id_idx = -1
        label_asym_id_idx = label_seq_id_idx = cartn_x_idx = cartn_y_idx = cartn_z_idx = -1

        for line in lines:
            line = line.strip()
            if line.startswith('loop_'):
                headers = []
                in_atom_site_loop = False
            elif line.startswith('_atom_site.'):
                header = line.split('.', 1)[1]
                headers.append(header)
                # Map header fields to indices
                if header == 'id': id_idx = len(headers)-1
                elif header == 'label_atom_id': label_atom_id_idx = len(headers)-1
                elif header == 'label_comp_id': label_comp_id_idx = len(headers)-1
                elif header == 'label_asym_id': label_asym_id_idx = len(headers)-1
                elif header == 'label_seq_id': label_seq_id_idx = len(headers)-1
                elif header == 'Cartn_x': cartn_x_idx = len(headers)-1
                elif header == 'Cartn_y': cartn_y_idx = len(headers)-1
                elif header == 'Cartn_z': cartn_z_idx = len(headers)-1
            elif headers and not line.startswith('_'):
                in_atom_site_loop = True
            else:
                in_atom_site_loop = False

            if in_atom_site_loop and line and not line.startswith(('#', 'loop_')):
                parts = line.split()
                try:
                    atom = Atom(
                        number=int(parts[id_idx]),
                        residue=parts[label_seq_id_idx] if label_seq_id_idx != -1 else '?',
                        type=parts[label_atom_id_idx],
                        chain=parts[label_asym_id_idx],
                        x=float(parts[cartn_x_idx]),
                        y=float(parts[cartn_y_idx]),
                        z=float(parts[cartn_z_idx])
                    )
                    self.atoms.append(atom)
                except (IndexError, ValueError):
                    continue 

    # def predict_crosslinks(self, residues_1: str, residues_2: str, 
    #                       min_length: float, max_length: float) -> 'CrossLinkDataset':
    #     # Implementation not shown for brevity
    #     pass

    def __iter__(self):
        return iter(self.atoms)

class Export:
    def to_chimerax_pseudobonds(
        self, 
        pcd: ProteinChainDataset, 
        folder_path: str, 
        file_name: str,  
        protein_structure: str = None,
        min_distance: float = 0,
        max_distance: float = sys.maxsize,
        diameter: float = 0.2, 
        dashes: int = 1,
        color_valid_distance: str = "48cae4", #Sky Blue
        color_invalid_outsider: str = "#d62828", # Red
    ) -> None:

        os.makedirs(folder_path, exist_ok=True)
        xl_frequencies: Set[int] = set(self.xls_site_count.values())

        def _write_to_pb_file(buffer: str, filename: str):
            print(f'Writing {len(buffer)} crosslinks to {filename}')
            if buffer:
                file_path = os.path.join(folder_path, filename)
                with open(file_path, "w") as file:
                    file.write(f"; dashes = {dashes}\n; radius = {diameter}\n{buffer}")

        for xl_frequency in xl_frequencies:
            buffer_homotypical_xl = []
            buffer_heterotypical_intra_xl = []
            buffer_heterotypical_inter_xl = []

            for key, value in self.xls_site_count.items():
                if value != xl_frequency:
                    continue

                crosslink = self._validate_terminus_sites(key)
                chains = pcd[crosslink.protein_1]

                if crosslink.is_homotypical:
                    for c1 in chains:
                        for c2 in chains:
                            if c1 != c2:  # Ensure unique chain pairing
                                buffer_homotypical_xl.append(f"/{c1}:{crosslink.num_site_1}@CA\t/{c2}:{crosslink.num_site_2}@CA\t{color_valid_distance}")

                elif crosslink.is_interprotein:
                    chain1, chain2 = pcd[crosslink.protein_1], pcd[crosslink.protein_2]
                    for c1 in chain1:
                        for c2 in chain2:
                            buffer_heterotypical_inter_xl.append(f"/{c1}:{crosslink.num_site_1}@CA\t/{c2}:{crosslink.num_site_2}@CA\t{color_valid_distance}")

                else:  # Intraprotein heterotypical
                    for c1 in chains:
                        for c2 in chains:
                            buffer_heterotypical_intra_xl.append(f"/{c1}:{crosslink.num_site_1}@CA\t/{c2}:{crosslink.num_site_2}@CA\t{color_valid_distance}")

            # Writing to files
            _write_to_pb_file("\n".join(buffer_heterotypical_inter_xl), f"{file_name}_heterotypical_interprotein_xl_{xl_frequency}_rep.pb")
            _write_to_pb_file("\n".join(buffer_heterotypical_intra_xl), f"{file_name}_heterotypical_intraprotein_xl_{xl_frequency}_rep.pb")
            _write_to_pb_file("\n".join(buffer_homotypical_xl), f"{file_name}_homotypical_xl_{xl_frequency}_rep.pb")

        print(f"DB files saved to {folder_path}")

    def export_ppis_for_gephi(
        self, 
        pcd: ProteinChainDataset, 
        folder_path: str, 
        file_name: str,
        blank_counter: bool = False
    ) -> None:
        save_path = self._validate_gephi_format(folder_path, file_name)
        node_buffer = dict()
        edge_buffer = dict()
    
        for key, value in pcd:
            for i in value:
                node_buffer[i] = key
    
        for xl in self.xls_site_count:
            for chain1, id1 in node_buffer.items():
                for chain2, id2 in node_buffer.items():
                    if chain1 == chain2:
                        continue
                    xl = self._validate_terminus_sites(xl)
                    if xl.protein_1 == id1 and xl.protein_2 == id2:
                        pair = tuple(sorted((chain1, chain2)))
                        edge_buffer[pair] = edge_buffer.get(pair, 0) + 1
    
        if blank_counter:
            edge_buffer = {k: 1 for k in edge_buffer}
    
        self._create_gexf(save_path, node_buffer, edge_buffer)

    def export_aais_for_gephi(
        self, 
        crosslinks: CrossLinkDataset,
        pcd: ProteinChainDataset, 
        folder_path: str, 
        file_name: str,
        blank_counter: bool = False  
    ) -> None:
        crosslinks = copy.deepcopy(crosslinks)
        save_path = self._validate_gephi_format(folder_path, file_name)
        node_buffer = dict()
        edge_buffer = dict()
        size_buffer = dict()

        # Calculate min/max sites for each protein
        for key, _ in pcd:
            min_site = sys.maxsize
            max_site = 0

            for xl in self.xls_site_count:
                validated_xl = self._validate_terminus_sites(xl)
                if validated_xl.protein_1 == key:
                    min_site = min(min_site, validated_xl.num_site_1)
                    max_site = max(max_site, validated_xl.num_site_1)
                if validated_xl.protein_2 == key:
                    min_site = min(min_site, validated_xl.num_site_2)
                    max_site = max(max_site, validated_xl.num_site_2)

            size_buffer[key] = (min_site, max_site)

        for key, value in pcd:
            for chain in value:
                first = size_buffer[key][0]
                last = size_buffer[key][1] + 1
                for i in range(first, last):
                    node_id = f'{chain}{i}_{key}'
                    node_label = f'{chain}_{i}'
                    node_buffer[node_id] = node_label

                for i in range(first, last - 1):
                    node_id = f'{chain}{i}_{key}'
                    next_node = f'{chain}{i+1}_{key}'
                    edge_buffer[(node_id, next_node)] = 1

        # Map chains to their protein IDs
        chain_buffer = {chain: key for key, chains in pcd for chain in chains}

        # Process crosslinks between chains
        for xl in self.xls_site_count:
            validated_xl = self._validate_terminus_sites(xl)
            for chain1, id1 in chain_buffer.items():
                for chain2, id2 in chain_buffer.items():
                    if chain1 == chain2:
                        continue

                    if (validated_xl.protein_1 == id1 and validated_xl.protein_2 == id2):
                        first_id = f'{chain1}{validated_xl.num_site_1}_{id1}'
                        second_id = f'{chain2}{validated_xl.num_site_2}_{id2}'
                        pair = tuple(sorted((first_id, second_id)))
                        edge_buffer[pair] = edge_buffer.get(pair, 0) + 1

        self._create_gexf(save_path, node_buffer, edge_buffer)

    def _validate_terminus_sites(self, crosslink: 'CrossLink') -> 'CrossLink':
        buffer = copy.deepcopy(crosslink)
        if '{' in buffer.site_1:
            buffer.num_site_1 += 1 
        if '}' in buffer.site_1:
            buffer.num_site_1 -= 1 
        if '{' in buffer.site_2:
            buffer.num_site_2 += 1 
        if '}' in buffer.site_2:
            buffer.num_site_2 -= 1 
        return buffer

    def _validate_gephi_format(self, folder_name: str, file_name: str) -> str:
        file_extension = os.path.splitext(file_name)[1]
        lower_extension = file_extension.lower()
        if lower_extension != ".gexf":
            raise ValueError(f'ERROR! Wrong file extension in {file_name}. Only ".gexf" format is supported')
        os.makedirs(folder_name, exist_ok=True)
        return os.path.join(folder_name, file_name)

    def _create_gexf(self, file_path: str, node_buffer, edge_buffer) -> None:
        gexf = ET.Element('gexf', version='1.3', xmlns='http://www.gexf.net/1.3')
        graph = ET.SubElement(gexf, 'graph', mode='static', defaultedgetype='undirected')
    
        nodes = ET.SubElement(graph, 'nodes')
        for node_id, label in node_buffer.items():
            ET.SubElement(nodes, 'node', {'id': str(node_id), 'label': str(label)})
    
        edges = ET.SubElement(graph, 'edges')
        for edge_id, ((source, target), weight) in enumerate(edge_buffer.items()):
            ET.SubElement(edges, 'edge', {
                'id': str(edge_id),
                'source': str(source),
                'target': str(target),
                'weight': str(weight)
            })
    
        rough_string = ET.tostring(gexf, encoding='utf-8')
        reparsed = minidom.parseString(rough_string)
        xml_content = reparsed.toprettyxml(indent="  ")
        xml_content = xml_content.replace('<?xml version="1.0" ?>',
                            '<?xml version="1.0" encoding="UTF-8"?>')

        with open(file_path, 'w', encoding='utf-8') as f:
            f.write(xml_content)

