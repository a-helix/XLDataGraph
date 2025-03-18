# import pytest
# from xldg.circos import *
# from xldg.xl import *
# from xldg.utils import *
# from xldg.fasta import *


# Domain()
# Domain_Dataset()
#     filter_by_fasta
# Circos_Plot
#     save
#     set_xls_colors
# # def test_positive_list_specified_type_files_from_folder():
# #     files = PathUtil.list_specified_type_files_from_folder(CWD, '.fasta')
# #     assert len(files) == 3

# # def test_exception_generate_custom_list_with_int_ranges():
# #     with pytest.raises(ValueError, match = "ERROR! start value is greater than end value"):
# #         custom_list = DatasetUtil.generate_custom_list_with_int_ranges((3, 1), (7, 5), (11, 9))

# # Current Working Directory
# # CWD = os.path.join(os.getcwd(), "tests", "utils_test")

# # # Test Data Folder
# # TDF = os.path.join(os.getcwd(), "tests", "test_data")

# # # Read dmn files from the folder
# # domain_folder_path = PathUtil.list_specified_type_files_from_folder(CWD, '.dmn')
# # all_domains = Domain_Dataset(domain_folder_path)
    
# # # Read fasta files from the folder
# # fasta_folder_path = PathUtil.list_specified_type_files_from_folder(CWD, '.fasta')
# # all_fastas = FastaDataset(fasta_folder_path, 'Uniprot')
    
# # # Read zhrm files from the folder
# # zhrm_folder_path = PathUtil.list_specified_type_files_from_folder(CWD, '.zhrm')   
# # sorted_zhrm_file_path = PathUtil.sort_filenames_by_first_integer(zhrm_folder_path)
# # folder_content =  DatasetUtil.read_merox_zhrm_files_from_path_list(sorted_zhrm_file_path, 'DSBU')

# # filtered_score_folder_content =  DatasetUtil.filter_all_by_score(folder_content, threshold=0)
    
# # combined_replicas = DatasetUtil.combine_replicas_in_CrossLinkDataset(filtered_score_folder_content, n=4)

# # dataset_1 = combined_replicas[0]
# # dataset_2 = combined_replicas[1]
# # dataset_3 = combined_replicas[2]

# # ### Individual datasets
# # config = Circos_Config(all_fastas, all_domains)

# # config.title = 'dataset_1'
# # plot = Circos_Plot(dataset_1, config)
# # plot.save(os.path.join(CWD, r'dataset_1.svg'))

# # config.title = 'dataset_2'
# # plot = Circos_Plot(dataset_2, config)
# # plot.save(os.path.join(CWD, r'dataset_2.svg'))

# # config.title = 'dataset_3'
# # plot = Circos_Plot(dataset_3, config)
# # plot.save(os.path.join(CWD, r'dataset_3.svg'))