import os
from xldg.utils import PathUtil

# Test cases for PathUtil
cwd = os.path.join(os.getcwd(), "tests", "utils_test")

def test_positive_list_specified_type_files_from_folder():
    files = PathUtil.list_specified_type_files_from_folder(cwd, '.fasta')
    assert len(files) == 3

def test_negative_list_specified_type_files_from_folder():
    files = PathUtil.list_specified_type_files_from_folder(cwd, '.none')
    assert len(files) == 0

def test_positive_sort_filenames_by_first_integer():
    files = PathUtil.list_specified_type_files_from_folder(cwd, '.fasta')
    sorted_files = PathUtil.sort_filenames_by_first_integer(files)
    reference = [os.path.join(cwd, '1_BSA.fasta'), 
                 os.path.join(cwd, '2_BSA.fasta'), 
                 os.path.join(cwd, '3_BSA.fasta')]
    assert sorted_files == reference

def test_negative_sort_filenames_by_first_integer():
    files = PathUtil.list_specified_type_files_from_folder(cwd, '.none')
    sorted_files = PathUtil.sort_filenames_by_first_integer(files)
    reference = []
    assert sorted_files == reference
