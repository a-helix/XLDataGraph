import pytest
import os
import random
from xldg.utils import PathUtil, DatasetUtil

### Test cases for PathUtil
# Current Working Directory
CWD = os.path.join(os.getcwd(), "tests", "test_data", "utils_test")

def test_positive_list_specified_type_files_from_folder():
    files = PathUtil.list_specified_type_files_from_folder(CWD, '.fasta')
    assert len(files) == 3

def test_negative_list_specified_type_files_from_folder():
    files = PathUtil.list_specified_type_files_from_folder(CWD, '.none')
    assert len(files) == 0

def test_exception_in_list_specified_type_files_from_folder():
    folder = "non_existent_folder"
    with pytest.raises(FileNotFoundError):
        PathUtil.list_specified_type_files_from_folder(folder, ".zhrm")
    

def test_positive_sort_filenames_by_first_integer():
    files = PathUtil.list_specified_type_files_from_folder(CWD, '.fasta')
    random.shuffle(files)
    sorted_files = PathUtil.sort_filenames_by_first_integer(files)
    reference = [os.path.join(CWD, 'BSA_1.fasta'),
                 os.path.join(CWD, '2_BSA.fasta'),
                 os.path.join(CWD, 'BSA_3.fasta')]
    assert sorted_files == reference

def test_negative_sort_filenames_by_first_integer():
    files = PathUtil.list_specified_type_files_from_folder(CWD, '.none')
    sorted_files = PathUtil.sort_filenames_by_first_integer(files)
    reference = []
    assert sorted_files == reference

def test_ignore_argument_in_sort_filenames_by_first_integer():
    files = PathUtil.list_specified_type_files_from_folder(CWD, '.txt')
    sorted_files = PathUtil.sort_filenames_by_first_integer(files, ignore = 'abcd1234_')
    reference = [os.path.join(CWD, 'abcd1234_file1.txt'),
                 os.path.join(CWD, 'abcd1234_file2.txt'),
                 os.path.join(CWD, 'abcd1234_file3.txt')]
    assert sorted_files == reference

### Test cases for DatasetUtil
# Test Data Folder
TDF = os.path.join(os.getcwd(), "tests", "test_data", "zhrm")

def test_read_merox_zhrm_files_from_path_list():
    zhrm_folder_path = PathUtil.list_specified_type_files_from_folder(TDF, '.zhrm')
    folder_content = DatasetUtil.read_merox_zhrm_files_from_path_list(zhrm_folder_path, 'DSBU')
    assert len(folder_content) == 3

def test_combine_replicas_in_xl_dataset():
    zhrm_folder_path = PathUtil.list_specified_type_files_from_folder(TDF, '.zhrm')
    folder_content = DatasetUtil.read_merox_zhrm_files_from_path_list(zhrm_folder_path, 'DSBU')
    combined_replicas = DatasetUtil.combine_replicas_in_xl_dataset(folder_content, 3)
    assert len(combined_replicas) == 1

def test_exception_combine_replicas_in_xl_dataset():
    zhrm_folder_path = PathUtil.list_specified_type_files_from_folder(TDF, '.zhrm')
    folder_content = DatasetUtil.read_merox_zhrm_files_from_path_list(zhrm_folder_path, 'DSBU')
    with pytest.raises(Exception, match = "ERROR! dataset size 3 is not mutiple to n=4"):
        combined_replicas = DatasetUtil.combine_replicas_in_xl_dataset(folder_content, 4)

def test_xls_preservation_after_replicas_in_xl_dataset():
    zhrm_folder_path = PathUtil.list_specified_type_files_from_folder(TDF, '.zhrm')
    folder_content = DatasetUtil.read_merox_zhrm_files_from_path_list(zhrm_folder_path, 'DSBU')
    combined_replicas = DatasetUtil.combine_replicas_in_xl_dataset(folder_content, 3)
    sum_of_xls = sum([len(data) for data in combined_replicas])
    assert len(combined_replicas[0]) == sum_of_xls

def test_negative_filter_all_results_by_score():
    zhrm_folder_path = PathUtil.list_specified_type_files_from_folder(TDF, '.zhrm')
    folder_content = DatasetUtil.read_merox_zhrm_files_from_path_list(zhrm_folder_path, 'DSBU')
    combined_replicas = DatasetUtil.combine_replicas_in_xl_dataset(folder_content, 3)
    before_filtering = DatasetUtil.filter_all_results_by_score(combined_replicas, 0)
    after_filtering = DatasetUtil.filter_all_results_by_score(combined_replicas, 0)
    assert len(before_filtering[0]) == len(after_filtering[0])

def test_positive_min_argument_filter_all_results_by_score():
    zhrm_folder_path = PathUtil.list_specified_type_files_from_folder(TDF, '.zhrm')
    folder_content = DatasetUtil.read_merox_zhrm_files_from_path_list(zhrm_folder_path, 'DSBU')
    combined_replicas = DatasetUtil.combine_replicas_in_xl_dataset(folder_content, 3)
    filtered_data = DatasetUtil.filter_all_results_by_score(combined_replicas, 150)
    assert len(filtered_data[0]) == 25

def test_positive_max_argument_filter_all_results_by_score():
    zhrm_folder_path = PathUtil.list_specified_type_files_from_folder(TDF, '.zhrm')
    folder_content = DatasetUtil.read_merox_zhrm_files_from_path_list(zhrm_folder_path, 'DSBU')
    combined_replicas = DatasetUtil.combine_replicas_in_xl_dataset(folder_content, 3)
    filtered_data = DatasetUtil.filter_all_results_by_score(combined_replicas, max_score=1)
    assert len(filtered_data[0]) == 66

def test_positive_min_and_max_arguments_filter_all_results_by_score():
    zhrm_folder_path = PathUtil.list_specified_type_files_from_folder(TDF, '.zhrm')
    folder_content = DatasetUtil.read_merox_zhrm_files_from_path_list(zhrm_folder_path, 'DSBU')
    combined_replicas = DatasetUtil.combine_replicas_in_xl_dataset(folder_content, 3)
    filtered_data = DatasetUtil.filter_all_results_by_score(combined_replicas, 120, 150)
    assert len(filtered_data[0]) == 36

def test_exception_arguments_filter_all_results_by_score():
    zhrm_folder_path = PathUtil.list_specified_type_files_from_folder(TDF, '.zhrm')
    folder_content = DatasetUtil.read_merox_zhrm_files_from_path_list(zhrm_folder_path, 'DSBU')
    combined_replicas = DatasetUtil.combine_replicas_in_xl_dataset(folder_content, 3)
    with pytest.raises(ValueError, match = "ERROR! max_score is smaller than min_score"):
        filtered_data = DatasetUtil.filter_all_results_by_score(combined_replicas, 1, 0)
    
def test_combine_all_datasets():
    zhrm_folder_path = PathUtil.list_specified_type_files_from_folder(TDF, '.zhrm')
    folder_content = DatasetUtil.read_merox_zhrm_files_from_path_list(zhrm_folder_path, 'DSBU')
    sum_of_xls = sum([len(data) for data in folder_content])
    combined_datasets = DatasetUtil.combine_all_datasets(folder_content)
    assert len(combined_datasets) == sum_of_xls

def test_positive_generate_custom_list_with_int_ranges():
    custom_list = DatasetUtil.generate_custom_list_with_int_ranges((1, 3), (5, 7), (9, 11))
    assert custom_list == [1, 2, 3, 5, 6, 7, 9, 10, 11]

def test_exception_generate_custom_list_with_int_ranges():
    with pytest.raises(ValueError, match = "ERROR! start value is greater than end value"):
        custom_list = DatasetUtil.generate_custom_list_with_int_ranges((3, 1), (7, 5), (11, 9))

def test_combine_selected_datasets():
    zhrm_folder_path = PathUtil.list_specified_type_files_from_folder(TDF, '.zhrm')
    folder_content = DatasetUtil.read_merox_zhrm_files_from_path_list(zhrm_folder_path, 'DSBU')
    custom_list = DatasetUtil.generate_custom_list_with_int_ranges((1, 2))
    combined_dataset = DatasetUtil.combine_selected_datasets(folder_content, custom_list)
    reference_dataset = folder_content[1] + folder_content[2]
    assert len(combined_dataset) == len(reference_dataset)

def test_exception_combine_selected_datasets():
    zhrm_folder_path = PathUtil.list_specified_type_files_from_folder(TDF, '.zhrm')
    folder_content = DatasetUtil.read_merox_zhrm_files_from_path_list(zhrm_folder_path, 'DSBU')
    custom_list = DatasetUtil.generate_custom_list_with_int_ranges((1, 3))
    with pytest.raises(IndexError, match = "ERROR! index 3 out of given dataset_list range"):
        combined_dataset = DatasetUtil.combine_selected_datasets(folder_content, custom_list)
