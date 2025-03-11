import collections
import colorsys
import math
from typing import List, Tuple, Iterator, Dict
import zipfile
import os
import re
import io
import copy


class Protein_Chain_ID_Dataset:
    def __init__(self, pcid_file_path: str):
        self.path = pcid_file_path
        self.pcids = {}
        self._assign_pcids(pcid_file_path)

    def _assign_pcids(self, path_to_pcid_file: str) -> None:
        try:
            with open(path_to_pcid_file, 'r') as file:
                for line in file:
                    splited_line = line.replace('\n', '').split(',')
                    self.pcids[splited_line[0]] = splited_line[1:]

        except FileNotFoundError:
            raise ValueError(f'Protein_Chain_ID_Dataset error: File at {path_to_pcid_file} was not found.')
        except Exception as e:
            raise ValueError(f'Protein_Chain_ID_Dataset error: {e}')

    def __len__(self):
        return len(self.pcids)

    def __iter__(self):
        return iter(self.pcids.items())

    def __getitem__(self, key):
        return self.pcids[key]

    def __next__(self):
        return next(iter(self.pcids.items()))

class XL:
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


class XL_Dataset:
    def __init__(self, xls: List['XL']):        
        self.xls = xls
        self._remove_decoy_xls()

        self.size = len(self.xls)
        self.xls_site_count = self._quantify_elements(xls)

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
        if not isinstance(other, XL_Dataset):
            return NotImplemented  # Return if other is not a XL_Dataset

        combined_xls = self.xls + other.xls
        combined_xls_site_count = self.xls_site_count

        for site, count in other.xls_site_count.items():
            if site not in combined_xls_site_count:
                combined_xls_site_count[site] = count 
            else:
                combined_xls_site_count[site] += count

        final = XL_Dataset(combined_xls)
        final.xls_site_count = combined_xls_site_count

        return final

    def __getitem__(self, index):
        return self.xls[index]

    def __len__(self):
        return self.size
    
    def filter_by_score(self, threshold: int):
        filtered_list = []
        
        for xl in self.xls:
            if xl.score >= threshold:
                filtered_list.append(xl)
        
        unique_filtered_list = set(filtered_list)
        filterered_xls_site_count = {}

        for xl, count in self.xls_site_count.items():
            if xl in unique_filtered_list:
                filterered_xls_site_count[xl] = count
                    
        self.xls = filtered_list
        self.size = len(self.xls)
        self.xls_site_count = filterered_xls_site_count

    def filter_by_min_xl_replica(self, min_xl_replica: int):
        filtered_xls_site_count = {}
        filtered_xls = []

        for xl1, count in self.xls_site_count.items():
            if count >= min_xl_replica:
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

    def _update_xls_data(self, xls: List['XL']) -> None:
        filtered_xls_site_count = {}
        for xl1 in xls:
            for xl2, count in self.xls_site_count.items():
                if xl1 == xl2:
                    filtered_xls_site_count[xl2] = count

        self.xls = xls
        self.size = len(self.xls)
        self.xls_site_count = filtered_xls_site_count

    def set_xls_site_count_to_one(self) -> None:
        for key in self.xls_site_count.keys():
            self.xls_site_count[key] = 1  
        
    def _quantify_elements(self, elements: List['XL']) -> Dict['XL', int]:
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
                if xl.protein_1.startswith('DEC_') or xl.protein_2.startswith('DEC_'):
                    continue
            buffer.append(xl)
        self.xls = buffer

    def export_xls_counters(self, path: str, file_name: str, separator: str = '\t'):
        file = os.path.join(path, file_name)
        os.makedirs(path, exist_ok=True)

        str_file = file
        header = f'protein_1{separator}peptide_1{separator}from_1{separator}to_1{separator}site_1{separator}protein_2{separator}peptide_2{separator}from_2{separator}to_2{separator}site_2{separator}interprotein{separator}homotypical{separator}replicas\n'

        with open(file, 'w') as file:
            file.write(header)
            for xl, frequency in self.xls_site_count.items():
                file.write(f'{xl.protein_1}{separator}{xl.peptide_1}{separator}{xl.from_1}{separator}{xl.to_1}{separator}{xl.site_1}{separator}{xl.protein_2}{separator}{xl.peptide_2}{separator}{xl.from_2}{separator}{xl.to_2}{separator}{xl.site_2}{separator}{xl.is_interprotein}{separator}{xl.is_homotypical}{separator}{frequency}\n')

    def export_for_chimerax(self, path: str, name: str, pcid: Protein_Chain_ID_Dataset, diameter: int = 0.2, color_heterotypical_intraprotein_xl: str = '#21a2ed', color_heterotypical_interprotein_xl: str = '#00008B', color_homotypical_xl: str = '#ed2b21') -> None:
        new_folder = os.path.join(path, name)
        os.makedirs(new_folder, exist_ok=True)
        
        xl_frequencies = set(self.xls_site_count.values())

        for xl_frequency in xl_frequencies:
            parameters = f'; dashes = 1\n; radius = {diameter}\n'
            buffer_heterotypical_INTRAprotein_xl = ''
            buffer_heterotypical_INTERprotein_xl = ''
            buffer_homotypical_xl = ''
            
            for key, value in self.xls_site_count.items():
                if value == xl_frequency:
                    if key.is_homotypical:
                        chains = pcid[key.protein_1]
                        for c1 in chains:
                            for c2 in chains:
                                # ChimeraX 1.8 doesn't render the whole file whithout this check
                                if c1 != c2: 
                                    buffer_homotypical_xl += f'/{c1}:{key.num_site_1}@CA\t/{c2}:{key.num_site_2}@CA\t{color_homotypical_xl}\n'
                               
                    if not key.is_homotypical:
                        if key.is_interprotein:
                            chain1 = pcid[key.protein_1]
                            chain2 = pcid[key.protein_2]

                            for c1 in chain1:
                                for c2 in chain2:
                                    buffer_heterotypical_INTERprotein_xl += f'/{c1}:{key.num_site_1}@CA\t/{c2}:{key.num_site_2}@CA\t{color_heterotypical_interprotein_xl}\n'
                        else:
                            chains = pcid[key.protein_1]

                            for c1 in chains:
                                for c2 in chains:
                                    buffer_heterotypical_INTRAprotein_xl += f'/{c1}:{key.num_site_1}@CA\t/{c2}:{key.num_site_2}@CA\t{color_heterotypical_intraprotein_xl}\n'


            file_path = f'{new_folder}\\{name}_heterotypical_interaprotein_xl_{str(xl_frequency)}.pb'

            if buffer_heterotypical_INTERprotein_xl != '':
                with open(file_path, 'w') as file:
                    buffer_heterotypical_INTERprotein_xl = parameters + buffer_heterotypical_INTERprotein_xl
                    file.write(buffer_heterotypical_INTERprotein_xl)

            file_path = f'{new_folder}\\{name}_heterotypical_intraprotein_xl_{str(xl_frequency)}.pb'

            if buffer_heterotypical_INTRAprotein_xl != '':
                with open(file_path, 'w') as file:
                    buffer_heterotypical_INTRAprotein_xl = parameters + buffer_heterotypical_INTRAprotein_xl
                    file.write(buffer_heterotypical_INTRAprotein_xl)
                
            file_path = f'{new_folder}\\{name}_homotypical_xl_{str(xl_frequency)}.pb'

            if buffer_homotypical_xl != '':
                with open(file_path, 'w') as file:
                    buffer_homotypical_xl = parameters + buffer_homotypical_xl
                    file.write(buffer_homotypical_xl)

        print(f'DB files saved to {new_folder}')

    def export_for_alphalink(self, folder_path: str, file_name: str, FDR: float = 0.05, min_xl_replica: int = 1) -> None:
        xl_sites = set(self.xls_site_count.keys())
        buffer = ""

        site_1 = 0
        site_2 = 0

        for xl, count in self.xls_site_count.items():
            if count < min_xl_replica:
                continue
            site_1 = xl.num_site_1
            if "{" in xl.site_1:
                site_1 += 1

            site_2 = xl.num_site_2
            if "}" in xl.site_2:
                site_2 -= 1

            buffer += f'{site_1}\t{site_2}\t{FDR}\n'

        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
        
        file_path = os.path.join(folder_path, file_name)
        with open(file_path, "w") as file:
            file.write(buffer)
        
        print(f'Alphalink crosslink file saved to {file_path}')

    @classmethod
    def _filter_xls_site_count_by_list_of_xls(cls, xls_site_count: Dict['XL', int], list_of_xls: List['XL']) -> 'XL_Dataset':
        filtered_count = {k: xls_site_count[k] for k in xls_site_count if k in list_of_xls}

        return filtered_count

    @classmethod
    def combine_datasets(cls, datasets: List['XL_Dataset']) -> 'XL_Dataset':
        combined_xls = None
        for dataset in datasets:
            if combined_xls is None:
                combined_xls = dataset
            else:
                combined_xls += dataset

        return combined_xls
    
    @classmethod
    def unique_elements(cls, dataset1: 'XL_Dataset', dataset2: 'XL_Dataset') -> Tuple['XL_Dataset', 'XL_Dataset']:
        count1 = dataset1.xls_site_count
        count2 = dataset1.xls_site_count
        
        set1 = set(dataset1.xls)
        set2 = set(dataset2.xls)
        
        unique_to_dataset1 = [xl for xl in dataset1.xls if xl not in set2]
        unique_to_dataset2 = [xl for xl in dataset2.xls if xl not in set1]
        
        # Create datasets from unique elements
        unique_dataset1 = cls(unique_to_dataset1)
        unique_dataset2 = cls(unique_to_dataset2)
        
        # Set the xls_site_count for the unique datasets
        unique_dataset1.xls_site_count = XL_Dataset._filter_xls_site_count_by_list_of_xls(count1, unique_to_dataset1)
        unique_dataset2.xls_site_count = XL_Dataset._filter_xls_site_count_by_list_of_xls(count2, unique_to_dataset2)
        
        return unique_dataset1, unique_dataset2
    
    @classmethod
    def common_elements(cls, dataset1: 'XL_Dataset', dataset2: 'XL_Dataset') -> Tuple['XL_Dataset', 'XL_Dataset']:
        count1 = dataset1.xls_site_count
        count2 = dataset2.xls_site_count
        
        common_elements = set(dataset1.xls) & set(dataset2.xls)
        
        common_list1 = [xl for xl in dataset1.xls if xl in common_elements]
        common_list2 = [xl for xl in dataset2.xls if xl in common_elements]
        
        # Create datasets from common elements
        common_dataset1 = cls(common_list1)
        common_dataset2 = cls(common_list2)
        
        # Set the xls_site_count for the common datasets
        common_dataset1.xls_site_count = XL_Dataset._filter_xls_site_count_by_list_of_xls(count1, common_elements)
        common_dataset2.xls_site_count = XL_Dataset._filter_xls_site_count_by_list_of_xls(count2, common_elements)
        
        return common_dataset1, common_dataset2


class Fasta_Entity:
    def __init__(self, header: str, sequence: str, fasta_format: str):
        self.raw_header = header.replace('(', '').replace(')', '')  # Merox also removes scopes
        self.raw_sequence = sequence

        if fasta_format == 'Uniprot':
            self.db_id, self.prot_gene = self._split_uniprot_fasta_header(header)
        elif fasta_format == 'Araport11':
            self.db_id, self.prot_gene = self._split_araport11_fasta_header(header)
        elif fasta_format == 'Custom':
            self.db_id, self.prot_gene = self.raw_header, self.raw_header
        else:
            raise ValueError(f'Unknown FASTA format: {fasta_format}')

        self.sequence = '{' + self.raw_sequence + '}' #MeroX format of sequemce with N-term and C-term as figure brackets
        self.seq_length = len(self.sequence)

    def _split_uniprot_fasta_header(self, header: str) -> Tuple[str, str]:
        header = header.strip()

        # Split by '|', extracting UniProt ID
        splited_header = header.split('|')
        db_id = splited_header[1]

        # Extract 'GN=...' substring
        prot_gene_match = re.search(r'(GN=[^\s]+)', header)  # Match full 'GN=...'
        prot_gene = prot_gene_match.group(1) if prot_gene_match else ''  # Extract matched value

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
                self.prot_gene == other.prot_gene and
                self.sequence == other.sequence)

    def __hash__(self):
        return hash((self.raw_header,  
                     self.db_id, 
                     self.prot_gene,
                     self.sequence))

    def __lt__(self, other):
        return self.db_id < other.db_id
    
    def __gt__(self, other):
        return self.db_id > other.db_id


class Fasta_Dataset:
    def __init__(self, fasta_files_paths_list: List[str], fasta_format: str):
        self.fasta_format = fasta_format
        self.entities = self._extract_all_fasta_content_from_folder(fasta_files_paths_list)
        self.size = len(self.entities)
        self._index = 0  # Initialize an index for iteration
        
    def _extract_all_fasta_content_from_folder(self, fasta_files_path_list: List[str]) -> List['Fasta_Entity']:
        db_entities = []
        
        for file_path in fasta_files_path_list:
            raw_fasta_content = ''
            # Read content of all files
            try:
                with open(file_path, 'r') as file:
                    for line in file:
                        raw_fasta_content += line
            except FileNotFoundError:
                print(f'Fasta_Dataset error: File at {file_path} was not found.')
            except Exception as e:
                print(f'Fasta_Dataset error: {e}')
            
            # Separate sequences
            splited_fasta_content = raw_fasta_content.split('>')
            for fasta in splited_fasta_content:
                splited_fasta = fasta.split('\n')
                
                if len(splited_fasta) > 1:
                    header = '>' + splited_fasta[0]
                    sequence = ''.join(line.strip() for line in splited_fasta[1:])
                    db_entities.append(Fasta_Entity(header, sequence, self.fasta_format))

        sorted_db_entities = sorted(db_entities) # Unifies sector plotting order on a final figure
        return sorted_db_entities
    
    def __len__(self):
        return self.size

    def __iter__(self) -> Iterator[Fasta_Entity]:
        self._index = 0  # Reset index for new iteration
        return self
    
    def __next__(self) -> Fasta_Entity:
        if self._index < self.size:
            entity = self.entities[self._index]
            self._index += 1
            return entity
        else:
            raise StopIteration
        
    def remove_entities_without_merox_xl(self, merox_xls: 'XL_Dataset') -> None:
        filtered_entities = set()
        
        for fasta in self.entities:
            for xl in merox_xls:
                if xl.protein_1 == fasta.raw_header or xl.protein_2 == fasta.raw_header:
                    filtered_entities.add(fasta)
                    break  # Exit the inner loop if a match is found

        self.entities = sorted(list(filtered_entities)) # Unifies sector plotting order on a final figure
        self.size = len(self.entities)

    def find_protein_name_by_header_string(self, header: str) -> str:
        for fasta in self.entities:
            if fasta.raw_header == header:
                return fasta.prot_gene

    def extract_proteins_fasta_enteties_in_XL_Dataset(self, folder_path: str, file_name: str, merox_data: 'XL_Dataset') -> None:
        proteins_in_XL_Dataset = set()

        for xl, _ in merox_data.xls_site_count.items():
            proteins_in_XL_Dataset.add(xl.protein_1)
            proteins_in_XL_Dataset.add(xl.protein_2)

        print(len(proteins_in_XL_Dataset))
        text_output = ''

        print(proteins_in_XL_Dataset)

        for fasta in self.entities:
            print(f'|{fasta.raw_header}|')
            if fasta.raw_header in proteins_in_XL_Dataset:
                text_output += f'{fasta.raw_header}\n'
                text_output += f'{fasta.raw_sequence}\n'
        
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
        
        file_path = os.path.join(folder_path, file_name)
        with open(file_path, "w") as file:
            file.write(text_output)
        
        print(f'Fasta saved to {file_path}')


# TODO: Add class utilities
def list_specified_type_files_in_folder(folder_path: str, file_format: str) -> List[str]:
    files = []
    for file in os.listdir(folder_path):
        if file.endswith(file_format):
            files.append(os.path.join(folder_path, file))

    return files


def sort_by_first_integers_in_filename(strings: List[str]) -> List[str]:
    def extract_leading_integer_from_file_name(s: str) -> int:
        file_path_only = os.path.basename(s)
        match = re.match(r'^(\d+)_', file_path_only)
        return int(match.group(1)) if match else float('inf')

    sorted_strings = sorted(strings, key=extract_leading_integer_from_file_name)

    return sorted_strings


def _extract_merox_result_from_zhrm_file(path: str, linker: str) -> 'XL_Dataset':
    xls = []
    software = 'MeroX'

    with zipfile.ZipFile(path, 'r') as zip_ref:
        with zip_ref.open('Result.csv') as csv_file:
            for line in io.TextIOWrapper(csv_file, encoding='utf-8'):
                row = line.strip().split(';')
                xl = XL(row[7], row[6], row[8], row[9], row[20], 
                        row[11], row[10], row[12], row[13], row[21],
                        row[0], software, linker)
                xls.append(xl)

    dataset = XL_Dataset(xls)
    dataset.set_xls_site_count_to_one()

    return dataset


def read_merox_zhrm_files_from_path_list(path_list: List[str], linker: str = None) -> List['XL_Dataset']:
    file_content = []
    
    for path in path_list:
        print(f'Extracting: {path}')
        file_content.append(_extract_merox_result_from_zhrm_file(path, linker))   

    return file_content


def filter_all_results_by_score(dataset: List['XL_Dataset'], threshold: int = 50) -> List['XL_Dataset']:
    for data in dataset:
        data.filter_by_score(threshold)
    return dataset
    

def combine_replicas_in_xl_dataset(dataset: List['XL_Dataset'], n=3) -> List['XL_Dataset']:
    combined_dataset = []
    buffer = []
    
    if ((len(dataset) % n) != 0):
        raise Exception(f'ERROR! dataset size {len(dataset)} is not mutiple to n={n}')
    
    for data in dataset:
        if (len(buffer) == n):
            combined_dataset.append(XL_Dataset.combine_datasets(buffer))
            buffer.clear()
        buffer.append(data)
        
    combined_dataset.append(XL_Dataset.combine_datasets(buffer))

    return combined_dataset


def fuse_list_of_xl_datsets(dataset_list: List['XL_Dataset']) -> 'XL_Dataset':
    return XL_Dataset([element for sublist in dataset_list for element in sublist])


def generate_custom_list_with_int_ranges(*diapason) -> List[int]:
    custom_list = []

    for pair in diapason:
        start = pair[0]
        end = pair[1]
    
        for i in range(start, end):
            custom_list.append(i)

    return custom_list


def get_selected_xls_from_dataset_list(dataset_list: List['XL_Dataset'], indexes: List[int]) -> 'XL_Dataset':
    buffer = []
    
    for x in indexes:
        buffer.append(dataset_list[x])

    return XL_Dataset.combine_datasets(buffer) 