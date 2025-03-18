import zipfile
import os
import io
import re
import sys

from typing import List, Tuple

from xldg.xl import CrossLink, CrossLinkDataset


class PathUtil:
    @staticmethod
    def list_specified_type_files_from_folder(folder_path: str, file_format: str) -> List[str]:
        """List all files with a specific format from a folder."""
        if not os.path.isdir(folder_path):
            raise FileNotFoundError(f'Folder not found: {folder_path}')

        files = []
        for file in os.listdir(folder_path):
            if file.endswith(file_format):
                files.append(os.path.join(folder_path, file))

        return files

    @staticmethod
    def sort_filenames_by_first_integer(strings: List[str], ignore: str = None) -> List[str]:
        """Sort filenames by the first integer appearing in the filename."""
        def _extract_integer(s: str) -> int:
            file_name = os.path.basename(s)
            if ignore is not None:
                file_name = file_name.replace(ignore, '')

            match = re.search(r'(\d+)', file_name)
            return int(match.group(1)) if match else float('inf')

        return sorted(strings, key=_extract_integer)

class DatasetUtil:
    @staticmethod
    def _extract_merox_result_from_zhrm_file(path: str, linker: str) -> 'CrossLinkDataset':
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
            dataset.set_xls_site_count_to_one()
            return dataset

    @staticmethod
    def read_merox_zhrm_files_from_path_list(path_list: List[str], linker: str = None) -> List['CrossLinkDataset']:
        file_content = []
    
        for path in path_list:
            if not os.path.isfile(path):
                raise FileNotFoundError(f'File not found: {path}')
            print(f'Extracting: {path}')
            file_content.append(DatasetUtil._extract_merox_result_from_zhrm_file(path, linker))   
        return file_content

    @staticmethod
    def filter_all_by_score(dataset: List['CrossLinkDataset'], min_score: int = 0, max_score: int = sys.maxsize) -> List['CrossLinkDataset']:
        if  max_score < min_score:
            raise ValueError('ERROR! max_score is smaller than min_score')

        for data in dataset:
            data.filter_by_score(min_score, max_score)
        return dataset    

    @staticmethod
    def combine_replicas_in_xl_datasets(dataset: List['CrossLinkDataset'], n = 3) -> List['CrossLinkDataset']:
        combined_dataset = []
        buffer = []
    
        if ((len(dataset) % n) != 0):
            raise Exception(f'ERROR! dataset size {len(dataset)} is not mutiple to n={n}')
    
        for data in dataset:
            if (len(buffer) == n):
                combined_dataset.append(CrossLinkDataset.combine_datasets(buffer))
                buffer.clear()
            buffer.append(data)
        
        combined_dataset.append(CrossLinkDataset.combine_datasets(buffer))

        return combined_dataset

    @staticmethod
    def combine_all_datasets(dataset_list: List['CrossLinkDataset']) -> 'CrossLinkDataset':
        """Combines multiple CrossLinkDataset instances into a single CrossLinkDataset."""
        combined_xls = []
        for dataset in dataset_list:
            combined_xls.extend(dataset.xls)
        return CrossLinkDataset(combined_xls)

    @staticmethod
    def generate_custom_list_with_int_ranges(*diapason: Tuple[int, int]) -> List[int]:
        """Generates a list of integers by including the full range from start to end (inclusive)."""
        custom_list = []

        for start, end in diapason:
            if start > end:
                raise ValueError("ERROR! start value is greater than end value")
            custom_list.extend(range(start, end + 1))

        return custom_list

    @staticmethod
    def combine_selected_datasets(dataset_list: List['CrossLinkDataset'], indexes: List[int]) -> 'CrossLinkDataset':
        biggest_index = max(indexes)
        if biggest_index >= len(dataset_list):
            raise IndexError(f"ERROR! index {biggest_index} out of given dataset_list range")

        buffer = []
        for x in indexes:
            buffer.append(dataset_list[x])

        return CrossLinkDataset.combine_datasets(buffer) 