from dataclasses import dataclass
from typing import List, Tuple, Dict, Set, Iterator, Union
import os
import sys
import re
import copy

import xml.etree.ElementTree as ET
from xml.dom import minidom


class ProteinChainDataset:
    def __init__(self, pcd_data: str):
        self.pcds = {}
        self._assign_pcds(pcd_data)

    def _assign_pcds(self, pcd_data: str) -> None:
        for line in pcd_data.strip().split("\n"):
            splited_line = line.strip().split(",")

            if len(splited_line) < 2:
                raise ValueError("Wrong data format in provided input")

            self.pcds[splited_line[0]] = splited_line[1:]


    def __len__(self):
        return len(self.pcds)

    def __iter__(self):
        return iter(self.pcds.items())

    def __getitem__(self, key):
        return self.pcds[key]

    def __next__(self):
        return next(iter(self.pcds.items()))

class CrossLinkEntity:
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
    def __init__(self, xls: List['CrossLinkEntity']):
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

    def remove_interprotein_crosslinks(self):
        filtered_xls = []
        for xl in self.xls:
            if xl.is_homotypical:
                filtered_xls.append(xl)
                continue
            if xl.is_interprotein is False:
                filtered_xls.append(xl)

        self._update_xls_data(filtered_xls)

    def remove_intraprotein_crosslinks(self):
        filtered_xls = []
        for xl in self.xls:
            if xl.is_homotypical:
                filtered_xls.append(xl)
                continue
            if xl.is_interprotein is True:
                filtered_xls.append(xl)

        self._update_xls_data(filtered_xls)

    def remove_homotypic_crosslinks(self):
        filtered_xls = []
        for xl in self.xls:
            if xl.is_homotypical is False:
                filtered_xls.append(xl)

        self._update_xls_data(filtered_xls)

    def _update_xls_data(self, xls: List['CrossLinkEntity']) -> None:
        filtered_xls_site_count = {}
        for xl1 in xls:
            for xl2, count in self.xls_site_count.items():
                if xl1 == xl2:
                    filtered_xls_site_count[xl2] = count

        self.xls = xls
        self.size = len(self.xls)
        self.xls_site_count = filtered_xls_site_count

    def blank_crosslink_counter(self) -> None:
        for key in self.xls_site_count.keys():
            self.xls_site_count[key] = 1  
        
    def _quantify_elements(self, elements: List['CrossLinkEntity']) -> Dict['CrossLinkEntity', int]:
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

    def save_crosslink_counter(self, folder_path: str, file_name: str, separator: str = '\t'):
        file = os.path.join(folder_path, file_name)
        os.makedirs(folder_path, exist_ok=True)

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

    def save_crosslinks(self, folder_path: str, file_name: str, separator: str = '\t'):
        file = os.path.join(folder_path, file_name)
        os.makedirs(folder_path, exist_ok=True)

        str_file = file
        header = (
            f"software{separator}protein_1{separator}peptide_1{separator}from_1{separator}"
            f"to_1{separator}site_1{separator}protein_2{separator}peptide_2{separator}"
            f"from_2{separator}to_2{separator}site_2{separator}interprotein{separator}"
            f"homotypical{separator}score\n"
        )

        with open(file, 'w') as file:
            file.write(header)
            for xl in self.xls():
                file.write(
                    f"{xl.software}{separator}{xl.protein_1}{separator}{xl.peptide_1}{separator}"
                    f"{xl.from_1}{separator}{xl.to_1}{separator}{xl.site_1}{separator}"
                    f"{xl.protein_2}{separator}{xl.peptide_2}{separator}{xl.from_2}{separator}"
                    f"{xl.to_2}{separator}{xl.site_2}{separator}{xl.is_interprotein}{separator}"
                    f"{xl.is_homotypical}{separator}{xl.score}\n"
                )

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
            print(f'Writing {len(buffer)} CrossLinkEntitys to {filename}')
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

                CrossLinkEntity = self._validate_terminus_sites(key)
                chains = pcd[CrossLinkEntity.protein_1]

                if CrossLinkEntity.is_homotypical:
                    for c1 in chains:
                        for c2 in chains:
                            if c1 != c2:  # Ensure unique chain pairing
                                buffer_homotypical_xl.append(f"/{c1}:{CrossLinkEntity.num_site_1}@CA\t/{c2}:{CrossLinkEntity.num_site_2}@CA\t{color_valid_distance}")

                elif CrossLinkEntity.is_interprotein:
                    chain1, chain2 = pcd[CrossLinkEntity.protein_1], pcd[CrossLinkEntity.protein_2]
                    for c1 in chain1:
                        for c2 in chain2:
                            buffer_heterotypical_inter_xl.append(f"/{c1}:{CrossLinkEntity.num_site_1}@CA\t/{c2}:{CrossLinkEntity.num_site_2}@CA\t{color_valid_distance}")

                else:  # Intraprotein heterotypical
                    for c1 in chains:
                        for c2 in chains:
                            buffer_heterotypical_intra_xl.append(f"/{c1}:{CrossLinkEntity.num_site_1}@CA\t/{c2}:{CrossLinkEntity.num_site_2}@CA\t{color_valid_distance}")

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
        pcd: ProteinChainDataset, 
        folder_path: str, 
        file_name: str,
        blank_counter: bool = False  
    ) -> None:
        CrossLinkEntitys = copy.deepcopy(CrossLinkEntitys)
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

        # Process CrossLinkEntitys between chains
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

    def _validate_terminus_sites(self, CrossLinkEntity: 'CrossLinkEntity') -> 'CrossLinkEntity':
        buffer = copy.deepcopy(CrossLinkEntity)
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

        def _filter_xls_site_count(xls_site_count: Dict['CrossLinkEntity', int], list_of_xls: List['CrossLinkEntity']) -> 'CrossLinkDataset':
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

class FastaEntity:
    def __init__(self, header: str, sequence: str, fasta_format: str, remove_parenthesis: bool = False):
        self.raw_header = header
        self.remove_parenthesis = remove_parenthesis

        if self.remove_parenthesis:
            self.raw_header = header.replace('(', '').replace(')', '')  # Merox also removes scopes

        self.raw_sequence = sequence

        try:
            if fasta_format == 'Uniprot':
                self.db_id, self.prot_gene = self._split_uniprot_fasta_header(header)
            elif fasta_format == 'Araport11':
                self.db_id, self.prot_gene = self._split_araport11_fasta_header(header)
            elif fasta_format == 'Custom':
                self.db_id, self.prot_gene = self.raw_header, self.raw_header
            else:
                raise Exception()
        except Exception as e:
                raise ValueError(f'ERROR! Wrong FASTA format: {fasta_format}')

        #MeroX format of sequence with N-term and C-term as figure brackets
        self.sequence = '{' + self.raw_sequence + '}' 
        self.seq_length = len(self.sequence)

    def _split_uniprot_fasta_header(self, header: str) -> Tuple[str, str]:
        header = header.strip()

        splited_header = header.split('|')
        db_id = splited_header[1]

        prot_gene_match = re.search(r'(GN=[^\s]+)', header)  # Match full 'GN=...'
        prot_gene = prot_gene_match.group(1) if prot_gene_match else ''

        return db_id, prot_gene.replace('GN=', '')

    def _split_araport11_fasta_header(self, header: str) -> Tuple[str, str]:
        splited_header = header.strip().split('|')
        araport11_id = splited_header[0].replace(' ', '').replace('>', '')

        prot_gene = splited_header[1].replace("Symbols: ", "")
        prot_gene = prot_gene.split()
        prot_gene = prot_gene[0].replace(',', '').replace(' ', '')
        return araport11_id, prot_gene

    def __eq__(self, other):
        return (self.raw_header == other.raw_header and  
                self.db_id == other.db_id and 
                self.prot_gene == other.prot_gene)

    def __hash__(self):
        return hash((self.raw_header,  
                     self.db_id, 
                     self.prot_gene))

    def __lt__(self, other):
        return self.db_id < other.db_id
    
    def __gt__(self, other):
        return self.db_id > other.db_id

    def __str__(self) -> str:
        return f'{self.raw_header}\n{self.raw_sequence}'


class FastaDataset:
    def __init__(self, raw_fasta_content: str, fasta_format: str, remove_parenthesis: bool = False):
        self.fasta_format = fasta_format
        self.remove_parenthesis = remove_parenthesis
        self.entities = self._parse_fasta_content(raw_fasta_content)
        self.size = len(self.entities)
        self._iter_index = 0 

    def _parse_fasta_content(self, raw_fasta_content: str) -> List['FastaEntity']:
        db_entities: Set['FastaEntity'] = set()

        if raw_fasta_content:
            db_entities.update(self._process_fasta_entries(raw_fasta_content))

        return sorted(db_entities)

    def _process_fasta_entries(self, raw_fasta_content: str) -> List['FastaEntity']:
        entries = []
        current_header = None
        current_sequence = []

        for line in raw_fasta_content.splitlines():
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_header:
                    entries.append(FastaEntity(current_header, "".join(current_sequence), self.fasta_format))
                current_header = line[1:].strip()
                current_sequence = []
            else:
                current_sequence.append(line)

        if current_header:
            entries.append(FastaEntity(current_header, "".join(current_sequence), self.fasta_format))

        return entries
    
    def __len__(self):
        return self.size

    def __iter__(self) -> Iterator['FastaEntity']:
        self._iter_index = 0
        return self
    
    def __next__(self) -> 'FastaEntity':
        if self._iter_index < self.size:
            entity = self.entities[self._iter_index]
            self._iter_index += 1
            return entity
        else:
            raise StopIteration

    def __str__(self) -> str:
        return "\n".join(str(fasta) for fasta in self.entities)
        
    def filter_by_crosslinks(self, merox_xls: 'CrossLinkDataset') -> None:
        filtered_entities = set()
        
        for fasta in self.entities:
            for xl in merox_xls:
                if xl.protein_1 == fasta.raw_header or xl.protein_2 == fasta.raw_header:
                    filtered_entities.add(fasta)
                    break

        # Unifies sector plotting order on a final figure
        self.entities = sorted(list(filtered_entities)) 
        self.size = len(self.entities)

    def find_gene_by_fasta_header(self, header: str) -> str:
        for fasta in self.entities:
            if fasta.raw_header == header:
                return fasta.prot_gene

    def save(self, folder_path: str, file_name: str) -> None:
        text_output = self.__str__()
        
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
        
        file_path = os.path.join(folder_path, file_name)
        with open(file_path, "w") as file:
            file.write(text_output)
        
        print(f'FastaEntity saved to {file_path}')


class DomainEntity:
    def __init__(self, input: str):
        self.splited_data = input.split(',')
        if len(self.splited_data) == 5:
            self.gene = self.splited_data[0].replace(' ', '')
            self.start = int(self.splited_data[1].replace(' ', ''))
            self.end = int(self.splited_data[2].replace(' ', ''))
            self.color = self.splited_data[3].replace(' ', '')
            self.name = self.splited_data[4].replace('\n', '')
            self.base_color = False
        elif len(self.splited_data) == 2:
            self.gene = self.splited_data[0].replace(' ', '')
            self.color = self.splited_data[1].replace(' ', '').replace('\n', '')
            self.base_color = True
        else:
            raise ValueError(f'Unknown domain format: {input}')


class DomainDataset:
    def __init__(self, domain_files_paths_list: List[str]):
        self.domains = self._extract_all_domain_content_from_folder(domain_files_paths_list)
        self._index = 0  # Initialize an index for iteration
        self._size = len(self.domains)

    def _extract_all_domain_content_from_folder(self, domain_files_paths_list: List[str]) -> List['DomainEntity']:
        domains = []
        
        for file_path in domain_files_paths_list:
            # Read content of all files
            try:
                with open(file_path, 'r') as file:
                    for line in file:
                        # Ignore comments and empty lines
                        if line[0] == '#' or not line.strip():
                            continue
                        domains.append(DomainEntity(line))
            except FileNotFoundError:
                print(f'Domain_Dataset error: File at {file_path} was not found.')
            except Exception as e:
                print(f'Domain_Dataset error: {e}')

        return domains

    def __len__(self):
        return len(self.domains)

    def __iter__(self) -> Iterator['DomainEntity']:
        self._index = 0  # Reset index for new iteration
        return self
    
    def __next__(self) -> 'DomainEntity':
        if self._index < self._size:
            domain = self.domains[self._index]
            self._index += 1
            return domain
        else:
            raise StopIteration

    def filter_by_fasta(self, FastaDataset: 'FastaDataset') -> None:
        filtered_domains = []
        for domain in self.domains:
            for fasta in FastaDataset:
                if domain.gene == fasta.prot_gene:
                    filtered_domains.append(domain)
                    break

        self.domains = filtered_domains


@dataclass
class Atom:
    number: int
    residue: str  # Residue identifier (number + insertion code if present)
    type: str     # (e.g., 'N', 'CA')
    chain: str    
    x: float
    y: float
    z: float

class ProteinStructureDataset:
    def __init__(self, file_path: str, format: str):
        self.file_path = file_path
        self.atoms: List[Atom] = []
        content = self._read_file()
        if format == '.pdb':
            self._parse_pdb(content)
        elif format == '.cif':
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

    # def predict_CrossLinkEntitys(self, residues_1: str, residues_2: str, 
    #                       min_length: float, max_length: float) -> 'CrossLinkDataset':
    #     # Implementation not shown for brevity
    #     pass

    def __iter__(self):
        return iter(self.atoms)




    

