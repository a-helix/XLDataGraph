# import pytest
# import os
# from xldg.utils import PathUtil, DatasetUtil

# ### Test cases for PathUtil

# ProteinChainDataset
# FileNotFoundError
# Exception

# CrossLinkDataset
# filter_by_score
# filter_by_min_xl_replica
# remove_interprotein_xls
# remove_intraprotein_xls
# remove_homotypic_xls
# set_xls_site_count_to_one
# export_xls_counters
# export_for_chimerax
# export_for_alphalink
# unique_elements
# common_elements
# combine_datasets
# # Current Working Directory
# # CWD = os.path.join(os.getcwd(), "tests", "test_data", "zhrm")
# # def test_negative_list_specified_type_files_from_folder():
# #     files = PathUtil.list_specified_type_files_from_folder(CWD, '.none')
# #     assert len(files) == 0

# # def test_exception_in_list_specified_type_files_from_folder():
# #     folder = "non_existent_folder"
# #     with pytest.raises(FileNotFoundError):
# #         PathUtil.list_specified_type_files_from_folder(folder, ".zhrm")