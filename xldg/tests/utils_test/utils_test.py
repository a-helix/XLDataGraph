import os
import random
from xldg.utils import PathUtil, DatasetUtil

### Test cases for PathUtil
# Current Working Directory
CWD = os.path.join(os.getcwd(), "tests", "utils_test")

def test_positive_list_specified_type_files_from_folder():
    files = PathUtil.list_specified_type_files_from_folder(CWD, '.fasta')
    assert len(files) == 3

def test_negative_list_specified_type_files_from_folder():
    files = PathUtil.list_specified_type_files_from_folder(CWD, '.none')
    assert len(files) == 0

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
TDF = os.path.join(os.getcwd(), "tests", "test_data")

def test_read_merox_zhrm_files_from_path_list():
    zhrm_folder_path = PathUtil.list_specified_type_files_from_folder(TDF, '.zhrm')
    sorted_zhrm_file_path = DatasetUtil.read_merox_zhrm_files_from_path_list(zhrm_folder_path, 'DSBU')
    assert len(sorted_zhrm_file_path) == 3

# DatasetUtil.read_merox_zhrm_files_from_path_list
# filter_all_results_by_score
# combine_replicas_in_xl_dataset
# fuse_list_of_xl_datsets
# generate_custom_list_with_int_ranges
# combine_selected_datasets
