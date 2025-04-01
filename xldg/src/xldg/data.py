import zipfile
import os
import io
import re
import sys
import copy
from typing import List, Tuple, Union, Optional

from xldg.core import FastaDataset, DomainDataset, CrossLinkEntity, CrossLinkDataset, ProteinChainDataset, ProteinStructureDataset


class Path:
    @staticmethod
    def list_given_type_files(folder_path: str, file_format: str) -> List[str]:
        files = []
        for file in os.listdir(folder_path):
            if file.endswith(file_format):
                files.append(os.path.join(folder_path, file))

        return files

    @staticmethod
    def sort_filenames_by_first_integer(strings: List[str], ignore: Optional[str] = None) -> List[str]:
        def _extract_integer(s: str) -> int:
            file_name = os.path.basename(s)
            if ignore is not None:
                file_name = file_name.replace(ignore, '')

            match = re.search(r'(\d+)', file_name)
            return int(match.group(1)) if match else float('inf')

        return sorted(strings, key=_extract_integer)

    @staticmethod
    def validate_file_existence(path: str) -> None:
        if not os.path.isfile(path):
                raise FileNotFoundError(f'File not found: {path}')

    @staticmethod
    def validate_folder_existence(path: str) -> None:
        if not os.path.isdir(path):
                raise FileNotFoundError(f'Folder not found: {path}')

    @staticmethod
    def confirm_file_format(path: str, *formats: str) -> str:
        _, extension = os.path.splitext(path)
        extension = extension.lower().lstrip('.')

        expected_formats = {fmt.lower() for fmt in formats}

        # Check if the extension is in the expected formats
        if extension not in expected_formats:
            raise ValueError(f"Invalid file format. Expected {', '.join(expected_formats).upper()}, got {extension.upper()}")

        return extension

    @staticmethod
    def read_to_string(path: str) -> str:
        string_buffer = '';
        with open(path, 'r') as file:
            for line in file:
                string_buffer += line
        return string_buffer

class Fasta:
    @staticmethod
    def load_data(path: Union[str, List[str]], fasta_format: str, remove_parenthesis: bool = False) -> 'FastaDataset':
        if isinstance(path, str):
            path = [path] 
        
        all_contents = [] 
        
        for file_path in path:
            Path.validate_file_existence(file_path)  
            Path.confirm_file_format(file_path, 'fasta', 'fas') 
            all_contents.append(Path.read_to_string(file_path)) 

        combined_content = "\n".join(all_contents)  

        return FastaDataset(combined_content, fasta_format, remove_parenthesis)

    @staticmethod
    def filter_by_crosslinks(fasta: FastaDataset, crosslinks: CrossLinkDataset) -> FastaDataset:
        fasta_copy = copy.deepcopy(fasta)
        return fasta_copy.filter_by_crosslinks(fasta, crosslinks)

class Domain:
    @staticmethod
    def load_data(path: Union[str, List[str]]) -> 'DomainDataset':
        if isinstance(path, str):
            path = [path] 
        
        all_contents = [] 
        
        for file_path in path:
            Path.validate_file_existence(file_path)  
            Path.confirm_file_format(file_path, 'dmn') 
            all_contents.append(Path.read_to_string(file_path)) 

        combined_content = "\n".join(all_contents)  
        return DomainDataset(combined_content)

class MeroX:
    @staticmethod
    def load_data(path: Union[str, List[str]], linker: Optional[str] = None) -> Union['CrossLinkDataset', List['CrossLinkDataset']]:
        def process_file(file_path: str) -> 'CrossLinkDataset':
            Path.validate_file_existence(file_path)
            Path.confirm_file_format(file_path, 'zhrm')
            xls = []
            software = 'MeroX'

            with zipfile.ZipFile(file_path, 'r') as zip_ref:
                with zip_ref.open('Result.csv') as csv_file:
                    for line in io.TextIOWrapper(csv_file, encoding='utf-8'):
                        row = line.strip().split(';')
                        xl = CrossLinkEntity(row[7], row[6], row[8], row[9], row[20], 
                                       row[11], row[10], row[12], row[13], row[21],
                                       row[0], software, linker)
                        xls.append(xl)

            dataset = CrossLinkDataset(xls)
            dataset.blank_replica_counter()
            return dataset

        if isinstance(path, list):
            return [process_file(file) for file in path]
        else:
            return process_file(path) 

class CrossLink:
    @staticmethod
    def filter_by_score(
        dataset: Union['CrossLinkDataset', List['CrossLinkDataset']], 
        min_score: int = 0, 
        max_score: int = sys.maxsize
    ) -> Union['CrossLinkDataset', List['CrossLinkDataset']]:

        if max_score < min_score:
            raise ValueError('max_score is smaller than min_score')
        dataset = copy.deepcopy(dataset)

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
    def filter_by_replica(
        dataset: Union['CrossLinkDataset', List['CrossLinkDataset']], 
        min_replica: int = 0, 
        max_replica: int = sys.maxsize
    ) -> Union['CrossLinkDataset', List['CrossLinkDataset']]:

        if max_replica < min_replica:
            raise ValueError('max_replica is smaller than min_replica')
        dataset = copy.deepcopy(dataset)

        if isinstance(dataset, CrossLinkDataset):
            dataset.filter_by_replica(min_replica, max_replica)
            return dataset

        if isinstance(dataset, list):
            for data in dataset:
                if not isinstance(data, CrossLinkDataset):
                    raise TypeError(f'Expected CrossLinkDataset, got {type(data)}')
                data.filter_by_replica(min_replica, max_replica)
            return dataset

        raise TypeError(f'Expected CrossLinkDataset or List[CrossLinkDataset], got {type(dataset)}')

    @staticmethod
    def blank_replica(
        dataset: Union['CrossLinkDataset', List['CrossLinkDataset']], 
    ) -> Union['CrossLinkDataset', List['CrossLinkDataset']]:

        dataset = copy.deepcopy(dataset)
        if isinstance(dataset, CrossLinkDataset):
            dataset.blank_replica_counter()
            return dataset
    
        if isinstance(dataset, list):
            for data in dataset:
                if not isinstance(data, CrossLinkDataset):
                    raise TypeError(f'Expected CrossLinkDataset, got {type(data)}')
                data.blank_replica_counter()
            return dataset

        raise TypeError(f'Expected CrossLinkDataset or List[CrossLinkDataset], got {type(dataset)}')

    @staticmethod
    def remove_interprotein(
        dataset: Union['CrossLinkDataset', List['CrossLinkDataset']], 
    ) -> Union['CrossLinkDataset', List['CrossLinkDataset']]:

        dataset = copy.deepcopy(dataset)
        if isinstance(dataset, CrossLinkDataset):
            dataset.remove_interprotein_crosslinks()
            return dataset
    
        if isinstance(dataset, list):
            for data in dataset:
                if not isinstance(data, CrossLinkDataset):
                    raise TypeError(f'Expected CrossLinkDataset, got {type(data)}')
                data.remove_interprotein_crosslinks()
            return dataset

        raise TypeError(f'Expected CrossLinkDataset or List[CrossLinkDataset], got {type(dataset)}')
     
    @staticmethod
    def remove_intraprotein(
        dataset: Union['CrossLinkDataset', List['CrossLinkDataset']], 
    ) -> Union['CrossLinkDataset', List['CrossLinkDataset']]:

        dataset = copy.deepcopy(dataset)
        if isinstance(dataset, CrossLinkDataset):
            dataset.remove_intraprotein_crosslinks()
            return dataset
    
        if isinstance(dataset, list):
            for data in dataset:
                if not isinstance(data, CrossLinkDataset):
                    raise TypeError(f'Expected CrossLinkDataset, got {type(data)}')
                data.remove_intraprotein_crosslinks()
            return dataset

        raise TypeError(f'Expected CrossLinkDataset or List[CrossLinkDataset], got {type(dataset)}')

    @staticmethod
    def remove_homotypic(
        dataset: Union['CrossLinkDataset', List['CrossLinkDataset']], 
    ) -> Union['CrossLinkDataset', List['CrossLinkDataset']]:

        dataset = copy.deepcopy(dataset)
        if isinstance(dataset, CrossLinkDataset):
            dataset.remove_homotypic_crosslinks()
            return dataset
    
        if isinstance(dataset, list):
            for data in dataset:
                if not isinstance(data, CrossLinkDataset):
                    raise TypeError(f'Expected CrossLinkDataset, got {type(data)}')
                data.remove_homotypic_crosslinks()
            return dataset

        raise TypeError(f'Expected CrossLinkDataset or List[CrossLinkDataset], got {type(dataset)}')

    @staticmethod
    def combine_replicas(dataset: List['CrossLinkDataset'], n = 3) -> List['CrossLinkDataset']:
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
    def combine_all(dataset_list: List['CrossLinkDataset']) -> 'CrossLinkDataset':
        combined_xls = None
        for dataset in dataset_list:
            if combined_xls is None:
                combined_xls = dataset
            else:
                combined_xls += dataset
        return combined_xls

    @staticmethod
    def combine_selected(dataset_list: List['CrossLinkDataset'], indexes: List[int]) -> 'CrossLinkDataset':
        biggest_index = max(indexes)
        smallest_index = min(indexes)
        last_dataset_index = len(dataset_list) - 1

        if biggest_index > len(dataset_list):
            raise IndexError(f"Index {biggest_index} out of given dataset_list range 0 to {last_dataset_index}")
        if smallest_index < 0:
            raise IndexError(f"Index {smallest_index} out of given dataset_list range 0 to {last_dataset_index}")

        buffer = None
        for x in indexes:
            if buffer is None:
                buffer = dataset_list[x]
            else:
                buffer += dataset_list[x]

        return buffer 

    @staticmethod
    def get_common(
        first_dataset: 'CrossLinkDataset', 
        second_dataset: 'CrossLinkDataset'
                   ) -> Tuple['CrossLinkDataset', 'CrossLinkDataset']:
        return CrossLinkDataset.common_elements(first_dataset, second_dataset)

    @staticmethod
    def get_unique(
        first_dataset: 'CrossLinkDataset', 
        second_dataset: 'CrossLinkDataset'
                   ) -> Tuple['CrossLinkDataset', 'CrossLinkDataset']:
        return CrossLinkDataset.unique_elements(first_dataset, second_dataset)


class ProteinChain:
    @staticmethod
    def load_data(path: str) -> ' ProteinChainDataset':
        Path.validate_file_existence(path)
        Path.confirm_file_format(path, 'pcd')
        content = Path.read_to_string(path)
        return ProteinChainDataset(content);


class ProteinStructure:
    @staticmethod
    def load_data(path: str) -> ' ProteinStructureDataset':
        Path.validate_file_existence(path)
        extension = Path.confirm_file_format(path, 'pdb', 'cif')
        content = Path.read_to_string(path)
        return ProteinStructureDataset(content, extension);


class Util:
    @staticmethod
    def generate_list_of_integers(*diapason: Tuple[int, int]) -> List[int]:
        custom_list = []

        for start, end in diapason:
            if start > end:
                raise ValueError(f"Start value {start} is greater than end value {end}")
            custom_list.extend(range(start, end + 1))
        return custom_list