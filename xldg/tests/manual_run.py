from src.xldg.circos import *
from src.xldg.xl import *

folder_path = r'D:\xldg_test'
# Read dmn files from the folder
domain_folder_path = list_specified_type_files_in_folder(folder_path, '.dmn')
all_domains = Domain_Dataset(domain_folder_path)
    
# Read fasta files from the folder
fasta_folder_path = list_specified_type_files_in_folder(folder_path, '.fasta')
all_fastas = Fasta_Dataset(fasta_folder_path, 'Araport11')
    
# Read zhrm files from the folder
zhrm_folder_path = list_specified_type_files_in_folder(folder_path, '.zhrm')   
sorted_zhrm_file_path = sort_by_first_integers_in_filename(zhrm_folder_path)
folder_content = read_merox_zhrm_files_from_path_list(sorted_zhrm_file_path, 'DSBU')

filtered_score_folder_content = filter_all_results_by_score(folder_content, threshold=0)
    
combined_replicas = combine_replicas_in_xl_dataset(filtered_score_folder_content, n=4)

first = combined_replicas[0]
second = combined_replicas[1]

### Individual datasets
config = Circos_Config(all_fastas, all_domains)
config.legend = 'Some text to see test functionality'
config.title = 'First'


plot = Circos_Plot(first, config)
plot.save(os.path.join(folder_path, r'first.svg'))

config.title = 'Second'
plot = Circos_Plot(second, config)
plot.save(os.path.join(folder_path, r'second.svg'))