from typing import List, Tuple, Dict, Set
import os
import sys
import re
import copy

import xml.etree.ElementTree as ET
from xml.dom import minidom


class ProteinChainDataset:
    def __init__(self, pcd_file_path: str):
        self.path = pcd_file_path
        self.pcds = {}
        self._assign_pcds(pcd_file_path)

    def _assign_pcds(self, path_to_pcd_file: str) -> None:
        try:
            with open(path_to_pcd_file, 'r') as file:
                for line in file:
                    splited_line = line.replace('\n', '').replace(' ', '').split(',')

                    if len(splited_line) < 2:
                        raise ValueError(
                            f'ERROR! Wrong data format in {path_to_pcd_file}')

                    self.pcds[splited_line[0]] = splited_line[1:]

        except FileNotFoundError:
            raise FileNotFoundError(f'ERROR! File at {path_to_pcd_file} was not found.')

    def __len__(self):
        return len(self.pcds)

    def __iter__(self):
        return iter(self.pcds.items())

    def __getitem__(self, key):
        return self.pcds[key]

    def __next__(self):
        return next(iter(self.pcds.items()))

class CrossLink:
    def _remove_text_in_brackets(self, s: str) -> str:
        pattern = r'\(.*?\)'
        cleaned_string = re.sub(pattern, '', s)
        return cleaned_string

    def _initialize_merox_xl(self):
        self.protein_1 = self._remove_text_in_brackets( self.protein_1).replace('  ', ' ') #Fix for a strange Merox assignment
        self.from_1 = int(self.from_1)
        self.to_1 = int(self.to_1)
        self.num_site_1 = self.from_1 + int(self.site_1[1:])
        self.protein_2 = self._remove_text_in_brackets( self.protein_2).replace('  ', ' ') #Fix for a strange Merox assignment
        self.from_2 = int(self.from_2)
        self.to_2 = int(self.to_2)
        self.num_site_2 =  self.from_2 + int(self.site_2[1:])
        self.score = int(self.score)

    def __init__(self, 
                 protein_1: str, 
                 peptide_1: str, 
                 from_1: str, 
                 to_1: str, 
                 site_1: str, 
                 protein_2: str, 
                 peptide_2: str, 
                 from_2: str, 
                 to_2: str, 
                 site_2: str, 
                 score: str, 
                 software: str, 
                 linker: str):

        self.software = software 
        self.linker = linker
        # First peptide data
        self.protein_1 = protein_1 # str
        self.peptide_1 = peptide_1 # str
        self.from_1 = from_1 # int
        self.to_1 = to_1 # int
        self.site_1 = site_1 # str
        self.num_site_1 = None # int
        # Second peptide data
        self.protein_2 = protein_2 # str
        self.peptide_2 = peptide_2 # str
        self.from_2 = from_2 # int
        self.to_2 = to_2 # int
        self.site_2 = site_2 # str
        self.num_site_2 =  None # int
        # Additional info
        self.score = score # int

        if self.software == 'MeroX':
             self._initialize_merox_xl()
        else:
            raise Exception(f'{self.software} is not supported.')

        self.str_info = f'{self.protein_1},{self.site_1},{self.protein_2},{self.site_2},{self.score}'
        self.is_interprotein = (self.protein_1 != self.protein_2)
        self.is_homotypical = (self.protein_1 == self.protein_2 and (self.num_site_1 == self.num_site_2 
                                                                     or self.peptide_1 == self.peptide_2))
    
    def __eq__(self, other):
        return (self.protein_1 == other.protein_1 and
                self.peptide_1 == other.peptide_1 and
                self.site_1 == other.site_1 and
                self.protein_2 == other.protein_2 and
                self.peptide_2 == other.peptide_2 and
                self.site_2 == other.site_2)  
    
    def __hash__(self):
        return hash((self.protein_1, 
                     self.peptide_1, 
                     self.from_1, 
                     self.to_1, 
                     self.site_1, 
                     self.protein_2, 
                     self.peptide_2, 
                     self.from_2, 
                     self.to_2, 
                     self.site_2))
    
    def __str__(self):
        return self.str_info


class CrossLinkDataset:
    def __init__(self, xls: List['CrossLink']):        
        self.xls = xls
        self._remove_decoy_xls()

        self.size = len(self.xls)
        self.xls_site_count = self._quantify_elements(self.xls)

    def __iter__(self):
        self._index = 0
        return self

    def __next__(self):
        if self._index < len(self.xls):
            result = self.xls[self._index]
            self._index += 1
            return result
        else:
            raise StopIteration

    def __add__(self, other):
        if not isinstance(other, CrossLinkDataset):
            return NotImplemented  # Return if other is not a CrossLinkDataset

        combined_xls = self.xls + other.xls
        combined_xls_site_count = self.xls_site_count

        for site, count in other.xls_site_count.items():
            if site not in combined_xls_site_count:
                combined_xls_site_count[site] = count 
            else:
                combined_xls_site_count[site] += count

        final = CrossLinkDataset(combined_xls)
        final.xls_site_count = combined_xls_site_count

        return final

    def __getitem__(self, index):
        return self.xls[index]

    def __len__(self):
        return self.size
    
    def filter_by_score(self, min_score: int = 0, max_score: int = sys.maxsize):
        filtered_list = []

        for xl in self.xls:
            if xl.score >= min_score and xl.score <= max_score:
                filtered_list.append(xl)
        
        unique_filtered_list = set(filtered_list)
        filterered_xls_site_count = {}

        for xl, count in self.xls_site_count.items():
            if xl in unique_filtered_list:
                filterered_xls_site_count[xl] = count
                    
        self.xls = filtered_list
        self.size = len(self.xls)
        self.xls_site_count = filterered_xls_site_count

    def filter_by_replica(self, min_rep: int = 1, max_rep: int = sys.maxsize):
        if min_rep > max_rep:
            raise ValueError('ERROR! max_rep is smaller than min_rep')

        filtered_xls_site_count = {}
        filtered_xls = []

        for xl1, count in self.xls_site_count.items():
            if count >= min_rep and count <= max_rep:
                filtered_xls_site_count[xl1] = count
                for xl2 in self.xls:
                    if xl2 == xl1:
                        filtered_xls.append(xl2)

        self.xls_site_count = filtered_xls_site_count
        self.xls = filtered_xls
        self.size = len(self.xls)

    def remove_interprotein_xls(self):
        filtered_xls = []
        for xl in self.xls:
            if xl.is_homotypical:
                filtered_xls.append(xl)
                continue
            if xl.is_interprotein is False:
                filtered_xls.append(xl)

        self._update_xls_data(filtered_xls)

    def remove_intraprotein_xls(self):
        filtered_xls = []
        for xl in self.xls:
            if xl.is_homotypical:
                filtered_xls.append(xl)
                continue
            if xl.is_interprotein is True:
                filtered_xls.append(xl)

        self._update_xls_data(filtered_xls)

    def remove_homotypic_xls(self):
        filtered_xls = []
        for xl in self.xls:
            if xl.is_homotypical is False:
                filtered_xls.append(xl)

        self._update_xls_data(filtered_xls)

    def _update_xls_data(self, xls: List['CrossLink']) -> None:
        filtered_xls_site_count = {}
        for xl1 in xls:
            for xl2, count in self.xls_site_count.items():
                if xl1 == xl2:
                    filtered_xls_site_count[xl2] = count

        self.xls = xls
        self.size = len(self.xls)
        self.xls_site_count = filtered_xls_site_count

    def set_xls_counter_to_one(self) -> None:
        for key in self.xls_site_count.keys():
            self.xls_site_count[key] = 1  
        
    def _quantify_elements(self, elements: List['CrossLink']) -> Dict['CrossLink', int]:
        element_counts = {}
        for element in elements:
            if element not in element_counts:
                element_counts[element] = 1 
            else:
                element_counts[element] += 1 

        return element_counts

    def _remove_decoy_xls(self) -> None:
        buffer = []
        for xl in self.xls:
            if xl.software == 'MeroX':
                # Ignore MeroX decoy matches
                if 'DEC_' in xl.protein_1 or 'DEC_' in xl.protein_2:
                    continue
            buffer.append(xl)
        self.xls = buffer

    def export_xls_counters(self, path: str, file_name: str, separator: str = '\t'):
        file = os.path.join(path, file_name)
        os.makedirs(path, exist_ok=True)

        str_file = file
        header = (
            f"software{separator}protein_1{separator}peptide_1{separator}from_1{separator}"
            f"to_1{separator}site_1{separator}protein_2{separator}peptide_2{separator}"
            f"from_2{separator}to_2{separator}site_2{separator}interprotein{separator}"
            f"homotypical{separator}replicas\n"
        )

        with open(file, 'w') as file:
            file.write(header)
            for xl, frequency in self.xls_site_count.items():
                file.write(
                    f"{xl.software}{separator}{xl.protein_1}{separator}{xl.peptide_1}{separator}"
                    f"{xl.from_1}{separator}{xl.to_1}{separator}{xl.site_1}{separator}"
                    f"{xl.protein_2}{separator}{xl.peptide_2}{separator}{xl.from_2}{separator}"
                    f"{xl.to_2}{separator}{xl.site_2}{separator}{xl.is_interprotein}{separator}"
                    f"{xl.is_homotypical}{separator}{frequency}\n"
                )

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

    def export_for_chimerax(
        self, 
        pcd: ProteinChainDataset, 
        folder_path: str, 
        file_name: str,  
        diameter: float = 0.2, 
        dashes: int = 1,
        color_heterotypical_intraprotein_xl: str = "#21a2ed", 
        color_heterotypical_interprotein_xl: str = "#00008B", 
        color_homotypical_xl: str = "#ed2b21"
    ) -> None: #TODO" ADD fit model
        """
        Exports cross-links for visualization in ChimeraX.

        Parameters:
            pcd (ProteinChainDataset): A dataset mapping proteins to their chain identifiers.
            folder_path (str): The directory where the output files will be saved.
            file_name (str, optional): Base name for the output files. Defaults to an empty string.
            diameter (float, optional): The radius of the cross-link visualization in ChimeraX. Default is 0.2.
            dashes (int, optional): The number of dashes used for rendering cross-links in ChimeraX. Default is 1.
            color_heterotypical_intraprotein_xl (str, optional): Color for intra-protein heterotypical cross-links. Default is "#21a2ed".
            color_heterotypical_interprotein_xl (str, optional): Color for inter-protein heterotypical cross-links. Default is "#00008B".
            color_homotypical_xl (str, optional): Color for homotypical cross-links. Default is "#ed2b21".

        Returns:
            None: Writes cross-link data to `.pb` files for ChimeraX visualization.

        The output files are saved in the specified folder and contain cross-links formatted for ChimeraX.
        The file names follow this format:
            - {file_name}_heterotypical_interprotein_xl_{frequency}_rep.pb
            - {file_name}_heterotypical_intraprotein_xl_{frequency}_rep.pb
            - {file_name}_homotypical_xl_{frequency}_rep.pb
        """

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
                                buffer_homotypical_xl.append(f"/{c1}:{crosslink.num_site_1}@CA\t/{c2}:{crosslink.num_site_2}@CA\t{color_homotypical_xl}")

                elif crosslink.is_interprotein:
                    chain1, chain2 = pcd[crosslink.protein_1], pcd[crosslink.protein_2]
                    for c1 in chain1:
                        for c2 in chain2:
                            buffer_heterotypical_inter_xl.append(f"/{c1}:{crosslink.num_site_1}@CA\t/{c2}:{crosslink.num_site_2}@CA\t{color_heterotypical_interprotein_xl}")

                else:  # Intraprotein heterotypical
                    for c1 in chains:
                        for c2 in chains:
                            buffer_heterotypical_intra_xl.append(f"/{c1}:{crosslink.num_site_1}@CA\t/{c2}:{crosslink.num_site_2}@CA\t{color_heterotypical_intraprotein_xl}")

            # Writing to files
            _write_to_pb_file("\n".join(buffer_heterotypical_inter_xl), f"{file_name}_heterotypical_interprotein_xl_{xl_frequency}_rep.pb")
            _write_to_pb_file("\n".join(buffer_heterotypical_intra_xl), f"{file_name}_heterotypical_intraprotein_xl_{xl_frequency}_rep.pb")
            _write_to_pb_file("\n".join(buffer_homotypical_xl), f"{file_name}_homotypical_xl_{xl_frequency}_rep.pb")

        print(f"DB files saved to {folder_path}")
   
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
        pcd: ProteinChainDataset, 
        folder_path: str, 
        file_name: str,
        blank_counter: bool = False  
    ) -> None:

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

        # Generate nodes and consecutive edges
        for key, value in pcd:
            for chain in value:
                first = size_buffer[key][0]
                last = size_buffer[key][1] + 1
                for i in range(first, last):
                    node_id = f'{chain}{i}_{key}'
                    node_label = f'{chain}{i}'
                    node_buffer[node_id] = node_label

                # Add edges between consecutive residues (i → i+1)
                for i in range(first, last - 1):  # ✅ Stop at last-1
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

    @classmethod
    def unique_elements(cls, dataset1: 'CrossLinkDataset', dataset2: 'CrossLinkDataset') -> Tuple['CrossLinkDataset', 'CrossLinkDataset']:
        count1 = dataset1.xls_site_count
        count2 = dataset1.xls_site_count
        
        set1 = set(dataset1.xls)
        set2 = set(dataset2.xls)
        
        unique_to_dataset1 = [xl for xl in dataset1.xls if xl not in set2]
        unique_to_dataset2 = [xl for xl in dataset2.xls if xl not in set1]
        
        # Create datasets from unique elements
        unique_dataset1 = cls(unique_to_dataset1)
        unique_dataset2 = cls(unique_to_dataset2)

        def _filter_xls_site_count(xls_site_count: Dict['CrossLink', int], list_of_xls: List['CrossLink']) -> 'CrossLinkDataset':
            return {k: xls_site_count[k] for k in xls_site_count if k in list_of_xls}

        # Set the xls_site_count for the unique datasets
        unique_dataset1.xls_site_count = _filter_xls_site_count(count1, unique_to_dataset1)
        unique_dataset2.xls_site_count = _filter_xls_site_count(count2, unique_to_dataset2)
        
        return unique_dataset1, unique_dataset2



    @classmethod
    def common_elements(cls, dataset1: 'CrossLinkDataset', dataset2: 'CrossLinkDataset') -> Tuple['CrossLinkDataset', 'CrossLinkDataset']:
        count1 = dataset1.xls_site_count
        count2 = dataset2.xls_site_count
        
        common_elements = set(dataset1.xls) & set(dataset2.xls)
        
        common_list1 = [xl for xl in dataset1.xls if xl in common_elements]
        common_list2 = [xl for xl in dataset2.xls if xl in common_elements]
        
        common_dataset1 = cls(common_list1)
        common_dataset2 = cls(common_list2)

        common_dataset1.xls_site_count = {k: count1[k] for k in count1 if k in common_elements}
        common_dataset2.xls_site_count = {k: count2[k] for k in count2 if k in common_elements}
        
        return common_dataset1, common_dataset2

    
    @classmethod
    def combine_datasets(cls, datasets: List['CrossLinkDataset']) -> 'CrossLinkDataset':
        combined_xls = None
        for dataset in datasets:
            if combined_xls is None:
                combined_xls = dataset
            else:
                combined_xls += dataset

        return combined_xls
