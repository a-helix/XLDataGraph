import pytest
import os
import random
from xldg.utils import PathUtil, DatasetUtil

class TestPathUtil:
    @pytest.fixture(autouse=True)
    def setup(self):
        # Current Working Directory
        self.CWD = os.path.join(os.getcwd(), "tests", "test_data", "utils_test")
        self.none_files = PathUtil.list_specified_type_files_from_folder(self.CWD, '.none')
        self.fasta_files = PathUtil.list_specified_type_files_from_folder(self.CWD, '.fasta')

    def test_positive_list_specified_type_files_from_folder(self):
        assert len(self.fasta_files) == 3

    def test_negative_list_specified_type_files_from_folder(self):
        assert len(self.none_files) == 0

    def test_exception_in_list_specified_type_files_from_folder(self):
        folder = "non_existent_folder"
        with pytest.raises(FileNotFoundError):
            PathUtil.list_specified_type_files_from_folder(folder, ".zhrm")
    

    def test_positive_sort_filenames_by_first_integer(self):
        random.shuffle(self.fasta_files)
        sorted_files = PathUtil.sort_filenames_by_first_integer(self.fasta_files)
        reference = [os.path.join(self.CWD, 'BSA_1.fasta'),
                     os.path.join(self.CWD, '2_BSA.fasta'),
                     os.path.join(self.CWD, 'BSA_3.fasta')]
        assert sorted_files == reference

    def test_negative_sort_filenames_by_first_integer(self):
        sorted_files = PathUtil.sort_filenames_by_first_integer(self.none_files)
        reference = []
        assert sorted_files == reference

    def test_ignore_argument_in_sort_filenames_by_first_integer(self):
        files = PathUtil.list_specified_type_files_from_folder(self.CWD, '.txt')
        sorted_files = PathUtil.sort_filenames_by_first_integer(files, ignore = 'abcd1234_')
        reference = [os.path.join(self.CWD, 'abcd1234_file1.txt'),
                     os.path.join(self.CWD, 'abcd1234_file2.txt'),
                     os.path.join(self.CWD, 'abcd1234_file3.txt')]
        assert sorted_files == reference

class TestDatasetUtil:
    @pytest.fixture(autouse=True)
    def setup(self):
        # Current Working Directory
        self.CWD = os.path.join(os.getcwd(), "tests", "test_data", "utils_test")
        # Test Data Folder
        self.TDF = os.path.join(os.getcwd(), "tests", "test_data", "zhrm") 

        zhrm_folder_path = PathUtil.list_specified_type_files_from_folder(self.TDF, '.zhrm')
        self.folder_content = DatasetUtil.read_merox_zhrm_files_from_path_list(zhrm_folder_path, 'DSBU')
        self.combined_replicas = DatasetUtil.combine_replicas_in_xl_datasets(self.folder_content, 3)

    def test_read_merox_zhrm_files_from_path_list(self):
        assert len(self.folder_content) == 3

    def test_combine_replicas_in_xl_datasets(self):
        
        assert len(self.combined_replicas) == 1

    def test_exception_combine_replicas_in_xl_datasets(self):
        with pytest.raises(Exception, match = "ERROR! dataset size 3 is not mutiple to n=4"):
            combined_replicas = DatasetUtil.combine_replicas_in_xl_datasets(self.folder_content, 4)

    def test_xls_preservation_after_replicas_in_CrossLinkDataset(self):
        sum_of_xls = sum([len(data) for data in self.combined_replicas])
        assert len(self.combined_replicas[0]) == sum_of_xls

    def test_negative_filter_all_by_score(self):
        after_filtering = DatasetUtil.filter_all_by_score(self.combined_replicas, 0)
        assert len(self.combined_replicas[0]) == len(after_filtering[0])

    def test_positive_min_argument_filter_all_by_score(self):
        filtered_data = DatasetUtil.filter_all_by_score(self.combined_replicas, 150)
        assert len(filtered_data[0]) == 25

    def test_positive_max_argument_filter_all_by_score(self):
        filtered_data = DatasetUtil.filter_all_by_score(self.combined_replicas, max_score=1)
        assert len(filtered_data[0]) == 66

    def test_positive_min_and_max_arguments_filter_all_by_score(self):
        filtered_data = DatasetUtil.filter_all_by_score(self.combined_replicas, 120, 150)
        assert len(filtered_data[0]) == 36

    def test_exception_arguments_filter_all_by_score(self):
        with pytest.raises(ValueError, match = "ERROR! max_score is smaller than min_score"):
            filtered_data = DatasetUtil.filter_all_by_score(self.combined_replicas, 1, 0)
    
    def test_combine_all_datasets(self):
        sum_of_xls = sum([len(data) for data in self.folder_content])
        combined_datasets = DatasetUtil.combine_all_datasets(self.folder_content)
        assert len(combined_datasets) == sum_of_xls

    def test_positive_generate_custom_list_with_int_ranges(self):
        custom_list = DatasetUtil.generate_custom_list_with_int_ranges((1, 3), (5, 7), (9, 11))
        assert custom_list == [1, 2, 3, 5, 6, 7, 9, 10, 11]

    def test_exception_generate_custom_list_with_int_ranges(self):
        with pytest.raises(ValueError, match = "ERROR! start value is greater than end value"):
            custom_list = DatasetUtil.generate_custom_list_with_int_ranges((3, 1), (7, 5), (11, 9))

    def test_combine_selected_datasets(self):
        custom_list = DatasetUtil.generate_custom_list_with_int_ranges((1, 2))
        combined_dataset = DatasetUtil.combine_selected_datasets(self.folder_content, custom_list)
        reference_dataset = self.folder_content[1] + self.folder_content[2]
        assert len(combined_dataset) == len(reference_dataset)

    def test_exception_combine_selected_datasets(self):
        custom_list = DatasetUtil.generate_custom_list_with_int_ranges((1, 3))
        with pytest.raises(IndexError, match = "ERROR! index 3 out of given dataset_list range"):
            combined_dataset = DatasetUtil.combine_selected_datasets(self.folder_content, custom_list)
