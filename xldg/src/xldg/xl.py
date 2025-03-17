from typing import List, Tuple, Dict
import os
import sys
import re


class ProteinChainDataset:
    def __init__(self, pcd_file_path: str):
        self.path = pcd_file_path
        self.pcids = {}
        self._assign_pcids(pcd_file_path)

    def _assign_pcids(self, path_to_pcd_file: str) -> None:
        try:
            with open(path_to_pcd_file, 'r') as file:
                for line in file:
                    splited_line = line.replace('\n', '').split(',')
                    self.pcids[splited_line[0]] = splited_line[1:]

        except FileNotFoundError:
            raise ValueError(f'ERROR! File at {path_to_pcd_file} was not found.')
        except Exception as e:
            raise ValueError(f'ERROR! {e}')

    def __len__(self):
        return len(self.pcids)

    def __iter__(self):
        return iter(self.pcids.items())

    def __getitem__(self, key):
        return self.pcids[key]

    def __next__(self):
        return next(iter(self.pcids.items()))

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
    
    def filter_by_score(self, min_score: int, max_score: int):
        filtered_list = []
        if max_score is None:
            max_score = sys.maxsize
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

    def _update_xls_data(self, xls: List['CrossLink']) -> None:
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

    def export_for_chimerax(self, path: str, name: str, pcd: ProteinChainDataset, diameter: int = 0.2, color_heterotypical_intraprotein_xl: str = '#21a2ed', color_heterotypical_interprotein_xl: str = '#00008B', color_homotypical_xl: str = '#ed2b21') -> None:
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
                        chains = pcd[key.protein_1]
                        for c1 in chains:
                            for c2 in chains:
                                # ChimeraX 1.8 doesn't render the whole file whithout this check
                                if c1 != c2: 
                                    buffer_homotypical_xl += f'/{c1}:{key.num_site_1}@CA\t/{c2}:{key.num_site_2}@CA\t{color_homotypical_xl}\n'
                               
                    if not key.is_homotypical:
                        if key.is_interprotein:
                            chain1 = pcd[key.protein_1]
                            chain2 = pcd[key.protein_2]

                            for c1 in chain1:
                                for c2 in chain2:
                                    buffer_heterotypical_INTERprotein_xl += f'/{c1}:{key.num_site_1}@CA\t/{c2}:{key.num_site_2}@CA\t{color_heterotypical_interprotein_xl}\n'
                        else:
                            chains = pcd[key.protein_1]

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

    # def export_for_alphalink(self, folder_path: str, file_name: str, FDR: float = 0.05, min_xl_replica: int = 1) -> None:
    #     xl_sites = set(self.xls_site_count.keys())
    #     buffer = ""

    #     site_1 = 0
    #     site_2 = 0

    #     for xl, count in self.xls_site_count.items():
    #         if count < min_xl_replica:
    #             continue
    #         site_1 = xl.num_site_1
    #         if "{" in xl.site_1:
    #             site_1 += 1

    #         site_2 = xl.num_site_2
    #         if "}" in xl.site_2:
    #             site_2 -= 1

    #         buffer += f'{site_1}\t{site_2}\t{FDR}\n'

    #     if not os.path.exists(folder_path):
    #         os.makedirs(folder_path)
        
    #     file_path = os.path.join(folder_path, file_name)
    #     with open(file_path, "w") as file:
    #         file.write(buffer)
        
    #     print(f'Alphalink crosslink file saved to {file_path}')

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
        
        # Create datasets from common elements
        common_dataset1 = cls(common_list1)
        common_dataset2 = cls(common_list2)
        
        # Set the xls_site_count for the common datasets
        common_dataset1.xls_site_count = CrossLinkDataset._filter_xls_site_count_by_list_of_xls(count1, common_elements)
        common_dataset2.xls_site_count = CrossLinkDataset._filter_xls_site_count_by_list_of_xls(count2, common_elements)
        
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
