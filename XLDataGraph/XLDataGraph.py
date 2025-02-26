# VERSION 0.2.1

import collections
import colorsys
import math
from typing import List, Tuple, Iterator, Dict
import zipfile
import os
import re
import io
import copy

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from matplotlib_venn import venn2
from pycirclize import Circos


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


class Merox_XL:
    def _remove_text_in_brackets(self, s: str) -> str:
        pattern = r'\(.*?\)'
        cleaned_string = re.sub(pattern, '', s)
        return cleaned_string

    def __init__(self, protein_1, peptide_1, from_1, to_1, site_1, protein_2, peptide_2, from_2, to_2, site_2, score):
        # First peptide data
        self.protein_1 = self._remove_text_in_brackets(protein_1).replace('  ', ' ') #Fix for a strange Merox assignment
        self.peptide_1 = peptide_1
        self.from_1 = int(from_1)
        self.to_1 = int(to_1)
        self.site_1 = site_1
        self.num_site_1 = self.from_1 + int(self.site_1[1:])
        # Second peptide data
        self.protein_2 = self._remove_text_in_brackets(protein_2).replace('  ', ' ') #Fix for a strange Merox assignment
        self.peptide_2 = peptide_2
        self.from_2 = int(from_2)
        self.to_2 = int(to_2)
        self.site_2 = site_2
        self.num_site_2 =  self.from_2 + int(self.site_2[1:])
        # Additional info
        self.score = int(score)
        self.str_info = f'{self.protein_1},{self.site_1},{self.protein_2},{self.site_2},{self.score}'

        self.is_homotypical = (self.protein_1 == self.protein_2 and (self.num_site_1 == self.num_site_2 
                                                                     or self.peptide_1 == self.peptide_2))
        self.is_interprotein = (self.protein_1 != self.protein_2)

    
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


class Merox_Dataset:
    def __init__(self, xls: List['Merox_XL']):
        self.xls = xls
        self._reomove_decoy_xls()

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
        if not isinstance(other, Merox_Dataset):
            return NotImplemented  # Return if other is not a Merox_Dataset

        combined_xls = self.xls + other.xls
        combined_xls_site_count = self.xls_site_count

        for site, count in other.xls_site_count.items():
            if site not in combined_xls_site_count:
                combined_xls_site_count[site] = count 
            else:
                combined_xls_site_count[site] += count

        final = Merox_Dataset(combined_xls)
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

    def _update_xls_data(self, xls: List['Merox_XL']) -> None:
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
        
    def _quantify_elements(self, elements: List['Merox_XL']) -> Dict['Merox_XL', int]:
        element_counts = {}
        for element in elements:
            if element not in element_counts:
                element_counts[element] = 1 
            else:
                element_counts[element] += 1 

        return element_counts

    def _reomove_decoy_xls(self) -> None:
        buffer = list()
        for xl in self.xls:
            # Ignore Merox decoy matched
            if(xl.protein_1[:4] == 'DEC_' or xl.protein_2[:4] == 'DEC_'):
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

    def export_for_chimerax(self, path: str, name: str, pcid: Protein_Chain_ID_Dataset, color_heterotypical_intraprotein_xl: str = '#21a2ed', color_heterotypical_interprotein_xl: str = '#00008B', color_homotypical_xl: str = '#ed2b21') -> None:
        new_folder = os.path.join(path, name)
        os.makedirs(new_folder, exist_ok=True)
        
        xl_frequencies = set(self.xls_site_count.values())

        for xl_frequency in xl_frequencies:
            diameter = 0.15

            parameters = f'; dashes = 1\n; radius = {(xl_frequency * diameter):.2f}\n'
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
    def _filter_xls_site_count_by_list_of_xls(cls, xls_site_count: Dict['Merox_XL', int], list_of_xls: List['Merox_XL']) -> 'Merox_Dataset':
        filtered_count = {k: xls_site_count[k] for k in xls_site_count if k in list_of_xls}

        return filtered_count

    @classmethod
    def combine_datasets(cls, datasets: List['Merox_Dataset']) -> 'Merox_Dataset':
        combined_xls = None
        for dataset in datasets:
            if combined_xls is None:
                combined_xls = dataset
            else:
                combined_xls += dataset

        return combined_xls
    
    @classmethod
    def unique_elements(cls, dataset1: 'Merox_Dataset', dataset2: 'Merox_Dataset') -> Tuple['Merox_Dataset', 'Merox_Dataset']:
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
        unique_dataset1.xls_site_count = Merox_Dataset._filter_xls_site_count_by_list_of_xls(count1, unique_to_dataset1)
        unique_dataset2.xls_site_count = Merox_Dataset._filter_xls_site_count_by_list_of_xls(count2, unique_to_dataset2)
        
        return unique_dataset1, unique_dataset2
    
    @classmethod
    def common_elements(cls, dataset1: 'Merox_Dataset', dataset2: 'Merox_Dataset') -> Tuple['Merox_Dataset', 'Merox_Dataset']:
        count1 = dataset1.xls_site_count
        count2 = dataset2.xls_site_count
        
        common_elements = set(dataset1.xls) & set(dataset2.xls)
        
        common_list1 = [xl for xl in dataset1.xls if xl in common_elements]
        common_list2 = [xl for xl in dataset2.xls if xl in common_elements]
        
        # Create datasets from common elements
        common_dataset1 = cls(common_list1)
        common_dataset2 = cls(common_list2)
        
        # Set the xls_site_count for the common datasets
        common_dataset1.xls_site_count = Merox_Dataset._filter_xls_site_count_by_list_of_xls(count1, common_elements)
        common_dataset2.xls_site_count = Merox_Dataset._filter_xls_site_count_by_list_of_xls(count2, common_elements)
        
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
        
    def remove_entities_without_merox_xl(self, merox_xls: 'Merox_Dataset') -> None:
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

    def extract_proteins_fasta_enteties_in_merox_dataset(self, folder_path: str, file_name: str, merox_data: 'Merox_Dataset') -> None:
        proteins_in_merox_dataset = set()

        for xl, _ in merox_data.xls_site_count.items():
            proteins_in_merox_dataset.add(xl.protein_1)
            proteins_in_merox_dataset.add(xl.protein_2)

        print(len(proteins_in_merox_dataset))
        text_output = ''

        print(proteins_in_merox_dataset)

        for fasta in self.entities:
            print(f'|{fasta.raw_header}|')
            if fasta.raw_header in proteins_in_merox_dataset:
                text_output += f'{fasta.raw_header}\n'
                text_output += f'{fasta.raw_sequence}\n'
        
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
        
        file_path = os.path.join(folder_path, file_name)
        with open(file_path, "w") as file:
            file.write(text_output)
        
        print(f'Fasta saved to {file_path}')


class Domain:
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


class Domain_Dataset:
    def __init__(self, domain_files_paths_list: List[str]):
        self.domains = self._extract_all_domain_content_from_folder(domain_files_paths_list)
        self.size = len(self.domains)
        self._index = 0  # Initialize an index for iteration
        
    def _extract_all_domain_content_from_folder(self, domain_files_paths_list: List[str]) -> List['Domain']:
        domains = []
        
        for file_path in domain_files_paths_list:
            # Read content of all files
            try:
                with open(file_path, 'r') as file:
                    for line in file:
                        if line[0] == '#':
                            continue
                        domains.append(Domain(line))
            except FileNotFoundError:
                print(f'Domain_Dataset error: File at {file_path} was not found.')
            except Exception as e:
                print(f'Domain_Dataset error: {e}')

        return domains

    def __len__(self):
        return self.size

    def __iter__(self) -> Iterator['Domain']:
        self._index = 0  # Reset index for new iteration
        return self
    
    def __next__(self) -> 'Domain':
        if self._index < self.size:
            domain = self.domains[self._index]
            self._index += 1
            return domain
        else:
            raise StopIteration

    def filter_by_fasta(self, fasta_dataset: 'Fasta_Dataset') -> None:
        filtered_domains = []
        for domain in self.domains:
            for fasta in fasta_dataset:
                if domain.gene == fasta.prot_gene:
                    filtered_domains.append(domain)
                    break

        self.domains = filtered_domains
        self.size = len(self.domains)


class Circos_Config:
        def __init__(self, 
                 # File input 
                 fasta: Fasta_Dataset, 
                 domains: Domain_Dataset = None,
                 # Text input 
                 legend: str = None, 
                 title: str = None, 
                 # Figure configs 
                 lable_interval: int = 20, 
                 space_between_sectors: int = 5,
                 domain_legend_distance: float = 1.15,
                 xl_legend_distance: float = 1.3,
                 xl_counter_distance: float = -0.15,
                 legend_distance: float = -0.15,
                 # Font configs 
                 title_font_size: int = 14,
                 lable_font_size: int = 14,
                 legend_font_size: int = 14,
                 prot_id_font_size: int = 14,
                 # Figure elements plotting configs 
                 plot_all_proteins: bool = False,
                 plot_protein_ids = True,
                 plot_xls_counter: bool = True,
                 plot_domain_legend: bool = True,
                 # XL configs 
                 min_xl_replica: int = 1,
                 plot_interprotein_xls: bool = True,
                 plot_intraprotein_xls: bool = True,
                 plot_homotypical_xls: bool = True):

            self.fasta = fasta
            self.domains = domains
            self.legend = legend
            self.title = title
            self.lable_interval = lable_interval
            self.space_between_sectors = space_between_sectors
            self.domain_legend_distance = domain_legend_distance
            self.xl_legend_distance = xl_legend_distance  
            self.xl_counter_distance = xl_counter_distance
            self.legend_distance = legend_distance
            self.title_font_size = title_font_size
            self.lable_font_size = lable_font_size
            self.legend_font_size = legend_font_size
            self.prot_id_font_size = prot_id_font_size
            self.plot_all_proteins = plot_all_proteins
            self.plot_protein_ids = plot_protein_ids
            self.plot_xls_counter = plot_xls_counter
            self.plot_domain_legend = plot_domain_legend
            self.min_xl_replica = min_xl_replica
            self.plot_interprotein_xls = plot_interprotein_xls
            self.plot_intraprotein_xls = plot_intraprotein_xls
            self.plot_homotypical_xls = plot_homotypical_xls


class Circos_Plot:  
    def __init__(self, xls: Merox_Dataset, config: Circos_Config):
        self.config = copy.deepcopy(config)
        self.xls = copy.deepcopy(xls)

        self.xls.filter_by_min_xl_replica(self.config.min_xl_replica)

        if self.config.plot_interprotein_xls is False:
            self.xls.remove_interprotein_xls()
        if self.config.plot_intraprotein_xls is False:
            self.xls.remove_intraprotein_xls()
        if self.config.plot_homotypical_xls is False:
            self.xls.remove_homotypic_xls()


        self.fasta = copy.deepcopy(config.fasta)
        if config.plot_all_proteins is False:
            self.fasta.remove_entities_without_merox_xl(self.xls)
        
        self.domains = None
        if self.config.domains is not None:
            self.domains = copy.deepcopy(self.config.domains)
            self.domains.filter_by_fasta(self.fasta)

        self.fig = None
        
        self.sectors = {prot.prot_gene: prot.seq_length for prot in self.fasta}
        self.prot_colors = self._assign_colors()
        self.circos = Circos(self.sectors, space=self.config.space_between_sectors)

        # XL colors
        self.heterotypic_intraprotein_xl_color = '#21a2ed' # Blue
        self.heterotypic_interprotein_xl_color ='#00008B' # Dark Blue
        self.homotypic_xl_color = '#ed2b21' # Red
        self.general_xl_color = '#7d8082' # Grey
     
    def save(self, path: str) -> None:
        if len(self.xls) == 0:
            print(f'WARNING: No crosslinks detected! Aborted save to {path}')
            return

        folder_path = os.path.dirname(path)
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
        
        self._plot_sectors()
        self._plot_xls()
        
        if (self.config.legend is not None):
            self._plot_user_legend()
            
        if (self.domains is not None and len(self.domains) != 0):
            self._plot_domains()
            
        self._plot_xl_legend()
        
        if self.config.plot_xls_counter is True:
            self._plot_xls_counter()
        
        if (self.config.title is not None):
            self._plot_title()

        self.fig.savefig(path)
        plt.close(self.fig)
        print(f'Circos plot saved to {path}')
    

    def _assign_colors(self) -> None:
        prot_colors = {}
        i = 0
        if self.domains is None:
            length = len(self.sectors)
            new_colors = self._generate_summer_colors(length)
            for prot in self.sectors:
                prot_colors[prot] = new_colors[i]
                i += 1
        else:
            for prot in self.sectors:
                prot_colors[prot] = '#C0C0C0'
                for domain in self.domains:
                    if domain.base_color is False:
                        continue
                    if prot == domain.gene:
                        prot_colors[prot] = domain.color
                        break        
                
        return prot_colors
    
    def _generate_summer_colors(self, num_colors: int) -> List[str]:
         summer_colors = []
         hue = 0.0  # Start at red (Hue 0)

         # Generate summer colors
         for _ in range(num_colors):
             lightness = 0.7  # High lightness for vibrant colors
             saturation = 0.8  # High saturation for bright colors
             r, g, b = colorsys.hls_to_rgb(hue, lightness, saturation)
             r = int(r * 255)
             g = int(g * 255)
             b = int(b * 255)
             summer_colors.append(f'#{r:02x}{g:02x}{b:02x}')

             # Increment hue for next color (golden ratio to avoid repeating colors)
             hue = (hue + 0.618033988749895) % 1.0

             # Skip hues to focus on typical summer colors (yellow, green, blue, pink)
             if 0.15 < hue < 0.3 or 0.55 < hue < 0.7:
                 hue = (hue + 0.2) % 1.0

         return summer_colors
    
    def _plot_sectors(self) -> None:
        for sector in self.circos.sectors:
            track = sector.add_track((92, 100))
            track.axis(fc = self.prot_colors[sector.name])
            if self.config.plot_protein_ids:
                sector.text(sector.name, color = '#3A3B3C', r = 110, size = self.config.prot_id_font_size) # Text lable

            if self.domains != None:
                for domain in self.domains:
                    if domain.gene != sector._name or domain.base_color is True:
                        continue
                    track2 = sector.add_track((92, 100))
                    track2.rect(domain.start, domain.end, fc=domain.color)
            
            track._start += 1 # Remove zero lable of the plot
            track.xticks_by_interval(self.config.lable_interval, label_size = self.config.lable_font_size) # Lable step
            track._start -= 1

    def _plot_xls(self) -> None:
        for xl, site_count in self.xls.xls_site_count.items():
            xl_color = self.heterotypic_intraprotein_xl_color
            plane = 2

            protein_1 = self.fasta.find_protein_name_by_header_string(xl.protein_1)
            protein_2 = self.fasta.find_protein_name_by_header_string(xl.protein_2)
            if protein_1 == None or protein_2 == None:
                continue

            if xl.is_homotypical:
                xl_color = self.homotypic_xl_color
                plane = 3
            elif xl.is_interprotein:
                xl_color = self.heterotypic_interprotein_xl_color
            
            self.circos.link((protein_1, xl.num_site_1, xl.num_site_1), (protein_2, xl.num_site_2, xl.num_site_2), ec=xl_color, zorder=plane, lw=site_count)
        
        self.fig = self.circos.plotfig()
    
    def _plot_xls_counter(self) -> None:
        total_xls_sites = 0
        site_counter = {}
        
        for xl, site_count in self.xls.xls_site_count.items():
            protein_1 = self.fasta.find_protein_name_by_header_string(xl.protein_1)
            protein_2 = self.fasta.find_protein_name_by_header_string(xl.protein_2)
            if protein_1 == None or protein_2 == None:
                continue
            
            total_xls_sites += 1
            
            if site_count in site_counter:
                site_counter[site_count] += 1
            else:
                site_counter[site_count] = 1
                
        sorted_site_counter = dict(sorted(site_counter.items()))
        
        if total_xls_sites > 0:
            text_lable = f'Total unique XLs: {total_xls_sites}\n'
            
            for key, value in sorted_site_counter.items():
                ending = ''
                if key > 1:
                    ending = 's'
                    
                text_lable += f'{key} replica{ending} unique XLs: {value}\n'
            
            self.fig.text(self.config.xl_counter_distance, 0.98, text_lable, fontsize=self.config.legend_font_size, va='top', ha='left')
      
    def _plot_user_legend(self) -> None:
        if self.config.legend != None:
            self.fig.text(self.config.legend_distance, 0.00, self.config.legend, va='bottom', ha='left', fontsize=self.config.legend_font_size)
           
    def _plot_domains(self) -> None:
        domains = [
            {'color': domain.color, 'label': domain.name}
            for domain in self.domains
            if domain.base_color is False
        ]
        legend_patches = []
        reference_buffer = []
        for item in domains:
            reference = item['color'] + item['label']
            if(reference in reference_buffer):
                continue
            
            check = item['label'].replace(' ', '')
            if(check != ''):
                legend_patches.append(mpatches.Patch(facecolor=item['color'], label=item['label'], linewidth=0.5, edgecolor='#3A3B3C'))
                reference_buffer.append(reference)
        
        if self.config.plot_domain_legend is True and len(legend_patches) != 0:
            self.fig.legend(handles=legend_patches, loc='lower right', bbox_to_anchor=(self.config.domain_legend_distance, 0), fontsize=self.config.legend_font_size)
    
    def _plot_xl_legend(self) -> None:
        most_frequent_xl = 0
        exhist_interprotein_xl = False
        exhist_intraprotein_xl = False
        exhist_homotypcal_xl = False

        for xl, site_count in self.xls.xls_site_count.items():
            if most_frequent_xl < site_count:
                most_frequent_xl = site_count

            if xl.is_homotypical:
                exhist_homotypcal_xl = True
            elif xl.is_interprotein:
                exhist_interprotein_xl = True
            else:
                exhist_intraprotein_xl = True
                
            
        if most_frequent_xl == 0:
            return

        legend_info = []
        if exhist_intraprotein_xl is True and self.config.plot_intraprotein_xls is True:
            legend_info.append({'label': 'Intraprotein unique XLs', 'color': self.heterotypic_intraprotein_xl_color, 'linewidth': 2})

        if exhist_interprotein_xl is True and self.config.plot_interprotein_xls is True:
            legend_info.append({'label': 'Interprotein unique XLs', 'color': self.heterotypic_interprotein_xl_color, 'linewidth': 2}) 

        if exhist_homotypcal_xl is True and self.config.plot_homotypical_xls is True:
            legend_info.append({'label': 'Homotypic unique XLs', 'color': self.homotypic_xl_color, 'linewidth': 2})

        if self.config.min_xl_replica == 1:
            legend_info.append({'label': '1-replica unique XLs', 'color': self.general_xl_color, 'linewidth': 1})
        
        if most_frequent_xl > 1:
            for i in range(2, most_frequent_xl+1):
                if i < self.config.min_xl_replica:
                    continue

                legend_info.append({'label': f'{i}-replicas unique XLs', 'color': self.general_xl_color, 'linewidth': i}) 
        
        legend_handles = [Line2D([0], [0], color=info['color'], linewidth=info['linewidth'], label=info['label']) for info in legend_info]
        self.fig.legend(handles=legend_handles, loc='upper right', bbox_to_anchor=(self.config.xl_legend_distance, 1), fontsize=self.config.legend_font_size)
    
    def _plot_title(self) -> None:
        if self.config.title is not None:    
            self.fig.text(0.5, 1.05, self.config.title, ha='center', va='center', fontsize=self.title_font_size)
    
    def set_xls_colors(self, 
                      heterotypic_intraprotein_xl_color = '#21a2ed', 
                      heterotypic_interprotein_xl_color = '#00008B', 
                      homotypic_xl_color = '#ed2b21', 
                      general_xl_color = '#7d8082') -> None:

        self.heterotypic_intraprotein_xl_color = heterotypic_intraprotein_xl_color
        self.heterotypic_interprotein_xl_color = heterotypic_interprotein_xl_color
        self.homotypic_xl_color = homotypic_xl_color
        self.general_xl_color = general_xl_color

   
def list_specified_type_files_in_folder(folder_path: str, file_format: str) -> List[str]:
    files = []
    for file in os.listdir(folder_path):
        if file.endswith(file_format):
            files.append(os.path.join(folder_path, file))

    return files


def sort_strings_by_first_integers_in_filename(strings: List[str]) -> List[str]:
    def extract_leading_integer_from_file_name(s: str) -> int:
        file_path_only = os.path.basename(s)
        match = re.match(r'^(\d+)_', file_path_only)
        return int(match.group(1)) if match else float('inf')

    sorted_strings = sorted(strings, key=extract_leading_integer_from_file_name)

    return sorted_strings


def extract_merox_result_from_zhrm_file(path: str) -> 'Merox_Dataset':
    xls = []
    with zipfile.ZipFile(path, 'r') as zip_ref:
        with zip_ref.open('Result.csv') as csv_file:
            for line in io.TextIOWrapper(csv_file, encoding='utf-8'):
                row = line.strip().split(';')
                xl = Merox_XL(row[7], row[6], row[8], row[9], row[20],
                        row[11], row[10], row[12], row[13], row[21],
                        row[0])
                xls.append(xl)

    dataset = Merox_Dataset(xls)
    dataset.set_xls_site_count_to_one()

    return dataset
  

def read_all_zhrm_files_from_list(path_list: List[str]) -> List['Merox_Dataset']:
    file_content = []
    
    for path in path_list:
        print(f'Extracting: {path}')
        file_content.append(extract_merox_result_from_zhrm_file(path))   

    return file_content


def filter_all_results_by_score(dataset: List['Merox_Dataset'], threshold: int = 50) -> List['Merox_Dataset']:
    for data in dataset:
        data.filter_by_score(threshold)
    return dataset
    
def combine_replicas_in_merox_results_dataset(dataset: List['Merox_Dataset'], n=3) -> List['Merox_Dataset']:
    combined_dataset = []
    buffer = []
    
    if ((len(dataset) % n) != 0):
        print(f'ERROR! dataset size {len(dataset)} is not mutiple to n={n}')
    
    for data in dataset:
        if (len(buffer) == n):
            combined_dataset.append(Merox_Dataset.combine_datasets(buffer))
            buffer.clear()
        buffer.append(data)
        
    combined_dataset.append(Merox_Dataset.combine_datasets(buffer))

    return combined_dataset


def generate_custom_list_with_int_ranges(*diapason) -> List[int]:
    custom_list = []

    for pair in diapason:
        start = pair[0]
        end = pair[1]
    
        for i in range(start, end):
            custom_list.append(i)

    return custom_list


def fuse_all_xl_datsets_in_one_dataset(dataset_list: List['Merox_Dataset']) -> 'Merox_Dataset':
    return Merox_Dataset([element for sublist in dataset_list for element in sublist])


def get_selected_xls_from_dataset_list(dataset_list: List['Merox_Dataset'], indexes: List[int]) -> 'Merox_Dataset':
    buffer = []
    
    for x in indexes:
        buffer.append(dataset_list[x])

    return Merox_Dataset.combine_datasets(buffer)  

# class Datasets_Plot:  
#     def __init__(self, xls: List['Merox_Dataset'], config: Datasets_Plot_Config):
#         self.config = copy.deepcopy(config)

def build_ven2_diagram_of_xl_sites(save_path: str, xls_list1: 'Merox_Dataset', xls_list2: 'Merox_Dataset', label1: str, label2: str, title: str = None, size: int = 16) -> None:
    plt.figure(figsize=(10, 10), dpi=600)
    
    set1 = set([str(sublist) for sublist in xls_list1])
    set2 = set([str(sublist) for sublist in xls_list2])
    venn = venn2([set1, set2], (label1, label2))

    # Customize colors
    venn.get_patch_by_id('10').set_color('#9AE66E') # pastel green
    venn.get_patch_by_id('01').set_color('#FAF278') # pastel yellow
    venn.get_patch_by_id('11').set_color('#87D5F8') # pastel blue

    # Label the regions with the number of elements
    for subset in ('10', '01', '11'):
        if venn.get_label_by_id(subset):
            venn.get_label_by_id(subset).set_text(f'{venn.get_label_by_id(subset).get_text()}')

    # Customize font size
    for text in venn.set_labels:
        text.set_fontsize(size)

    for text in venn.subset_labels:
        if text:  # Check if the subset label is not None
            text.set_fontsize(size)
    if title != None:
        plt.title(title).set_fontsize(18)

    plt.savefig(save_path)
    print(f'Venn2 diagram saved to {save_path}')
