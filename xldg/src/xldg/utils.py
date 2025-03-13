import zipfile
import os
import io
import re

from typing import List, Tuple

from xldg.xl import XL, XL_Dataset



class PathUtil:
    @staticmethod
    def list_specified_type_files_from_folder(folder_path: str, file_format: str) -> List[str]:
        """List all files with a specific format from a folder."""
        files = []
        for file in os.listdir(folder_path):
            if file.endswith(file_format):
                files.append(os.path.join(folder_path, file))

        return files

    @staticmethod
    def sort_filenames_by_first_integer(strings: List[str]) -> List[str]:
        """Sort filenames by the first integer appearing in the filename."""
        print('HERE')
        def extract_integer(s: str) -> int:
            file_name = os.path.basename(s)
            match = re.search(r'(\d+)', file_name)
            print(int(match.group(1)))
            return int(match.group(1)) if match else float('inf')
        print('HERE')
        return sorted(strings, key=extract_integer)

class DatasetUtil:
    @staticmethod
    def _extract_merox_result_from_zhrm_file(self, path: str, linker: str) -> 'XL_Dataset':
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

    @staticmethod
    def read_merox_zhrm_files_from_path_list(self, path_list: List[str], linker: str = None) -> List['XL_Dataset']:
        file_content = []
    
        for path in path_list:
            print(f'Extracting: {path}')
            file_content.append(self._extract_merox_result_from_zhrm_file(path, linker))   

        return file_content

    @staticmethod
    def filter_all_results_by_score(self, dataset: List['XL_Dataset'], threshold: int = 50) -> List['XL_Dataset']:
        for data in dataset:
            data.filter_by_score(threshold)
        return dataset    

    @staticmethod
    def combine_replicas_in_xl_dataset(self, dataset: List['XL_Dataset'], n=3) -> List['XL_Dataset']:
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

    @staticmethod
    def fuse_list_of_xl_datsets(self, dataset_list: List['XL_Dataset']) -> 'XL_Dataset':
        return XL_Dataset([element for sublist in dataset_list for element in sublist])

    @staticmethod
    def generate_custom_list_with_int_ranges(self, *diapason: Tuple[int, int]) -> List[int]:
        custom_list = []

        for pair in diapason:
            start = pair[0]
            end = pair[1]
    
            for i in range(start, end):
                custom_list.append(i)

        return custom_list

    @staticmethod
    def combine_selected_datasets(self, dataset_list: List['XL_Dataset'], indexes: List[int]) -> 'XL_Dataset':
        buffer = []
    
        for x in indexes:
            buffer.append(dataset_list[x])

        return XL_Dataset.combine_datasets(buffer) 