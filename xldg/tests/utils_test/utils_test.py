import os
import random
from xldg.utils import PathUtil, DatasetUtil

# Test cases for PathUtil
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

# Test cases for DatasetUtil