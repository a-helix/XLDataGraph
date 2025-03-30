# import pytest
# import os
# import copy

# from xldg.data import Path, Crosslink, FastaDataset
# from xldg.graphics import Circos, CircosConfig, DomainDataset, Domain




# class TestCircos:
#     @pytest.fixture(autouse=True)
#     def setup(self):
#         # Current Working Directory
#         self.CWD = os.path.join(os.getcwd(), "tests", "test_data", "circos_test")
#         fasta_single_path = os.path.join(os.getcwd(), "tests", "test_data", "fasta", "fasta_combined.fasta")
#         self.fasta_dataset = FastaDataset([fasta_single_path], "Custom")

#         # Domain Files Directory
#         DFD = os.path.join(os.getcwd(), "tests", "test_data", "dmn")
#         domain_path = Path.list_given_type_files(DFD, '.dmn')
#         self.domains = DomainDataset(domain_path)

#         merox_data = os.path.join(os.getcwd(), "tests", "test_data", "zhrm") 
#         zhrm_folder_path = Path.list_given_type_files(merox_data, '.zhrm')
#         folder_content = DatasetUtil.read_all_merox_files(zhrm_folder_path, 'DSBU')
#         self.combined_data = DatasetUtil.combine_all_datasets(folder_content)

#         self.config = CircosConfig(self.fasta_dataset)
#         self.circos = Circos(self.combined_data, self.config)

#     def _read_file(self, file_path: str, delete: bool = False):
#         content = None
#         with open(file_path, "r", encoding="utf-8") as f:
#             content = f.read()

#         if delete: 
#             os.remove(file_path)
#         return content

#     def test_positive_save_basic(self):
#         config = CircosConfig(self.fasta_dataset, 
#                               plot_domain_legend = False, 
#                               plot_protein_ids = False, 
#                               plot_counter = False,
#                               plot_xl_legend = False)
#         circos = Circos(self.combined_data, config)

#         save_path = os.path.join(self.CWD, "circos_basic.svg")
#         before_exhist = os.path.isfile(save_path)
#         circos.save(save_path)

#         after_exhist = os.path.isfile(save_path)
#         after_content = self._read_file(save_path, True)
#         reference_path = os.path.join(self.CWD, "circos_reference_basic.svg")
#         reference_content = self._read_file(reference_path)

#         assert (
#             before_exhist == False and 
#             after_exhist == True and
#             len(after_content) == len(reference_content)
#             )

#     def test_positive_save_cosmetics(self):
#         basic_config = CircosConfig(
#             self.fasta_dataset, 
#             self.domains,
#             plot_domain_legend=False
#             )

#         extended_config = CircosConfig(
#                                  fasta = self.fasta_dataset, 
#                                  domains = self.domains,
#                                  legend = 'This is legend\nSome important information\nMore important information', 
#                                  title = 'Long title of the circos plot', 
#                                  lable_interval = 35, 
#                                  space_between_sectors = 10,
#                                  domain_legend_distance = 1.2,
#                                  xl_legend_distance = 1.4,
#                                  xl_counter_distance = -0.1,
#                                  legend_distance = -0.1,
#                                  title_font_size = 20,
#                                  ruler_font_size = 12,
#                                  legend_font_size = 16,
#                                  prot_font_size = 18
#                                  )

#         basic_circos = Circos(self.combined_data, basic_config)
#         extended_circos = Circos(self.combined_data, extended_config)

#         basic_path = os.path.join(self.CWD, "circos_test_minimal_cosmetics.svg")
#         extended_path = os.path.join(self.CWD, "circos_test_extended_cosmetics.svg")

#         basic_circos.save(basic_path)
#         extended_circos.save(extended_path)

#         basic_content = self._read_file(basic_path, True)
#         extended_content = self._read_file(extended_path, True)

#         ref_basic_path = os.path.join(self.CWD, "circos_reference_minimal_cosmetics.svg")
#         ref_extended_path = os.path.join(self.CWD, "circos_reference_extended_cosmetics.svg")

#         ref_basic_content = self._read_file(ref_basic_path)
#         ref_extended_content = self._read_file(ref_extended_path)
#         assert (
#             len(basic_content) == len(ref_basic_content) and
#             len(extended_content) == len(ref_extended_content)
#             )

#     def test_positive_save_selected_xls(self):   
#         self.homotypic_config = CircosConfig(self.fasta_dataset, 
#                                              plot_interprotein_xls = False, 
#                                              plot_intraprotein_xls = False, 
#                                              plot_homotypical_xls = True)

#         self.inter_config = CircosConfig(self.fasta_dataset, 
#                                          plot_interprotein_xls = True, 
#                                          plot_intraprotein_xls = False, 
#                                          plot_homotypical_xls = False)

#         self.intra_config = CircosConfig(self.fasta_dataset, 
#                                          plot_interprotein_xls = False, 
#                                          plot_intraprotein_xls = True, 
#                                          plot_homotypical_xls = False)

#         self.homotypic_circos = Circos(self.combined_data, self.homotypic_config)
#         self.inter_circos = Circos(self.combined_data, self.inter_config)
#         self.intra_circos = Circos(self.combined_data, self.intra_config)

#         homotypic_path = os.path.join(self.CWD, "circos_test_homotypic.svg")
#         inter_path = os.path.join(self.CWD, "circos_test_inter.svg")
#         intra_path = os.path.join(self.CWD, "circos_test_intra.svg")

#         self.homotypic_circos.save(homotypic_path)
#         self.inter_circos.save(inter_path)
#         self.intra_circos.save(intra_path)

#         homotypic_content = self._read_file(homotypic_path, True)
#         inter_content = self._read_file(inter_path, True)
#         intra_content = self._read_file(intra_path, True)

#         ref_homotypic_path = os.path.join(self.CWD, "circos_reference_homotypic.svg")
#         ref_inter_path = os.path.join(self.CWD, "circos_reference_inter.svg")
#         ref_intra_path = os.path.join(self.CWD, "circos_reference_intra.svg")

#         ref_homotypic_content = self._read_file(ref_homotypic_path)
#         ref_inter_content = self._read_file(ref_inter_path)
#         ref_intra_content = self._read_file(ref_intra_path)

#         assert (
#             len(homotypic_content) == len(ref_homotypic_content) and
#             len(inter_content) == len(ref_inter_content) and
#             len(intra_content) == len(ref_intra_content)
#             )

#     def test_positive_save_min_max_rep(self):   
#         self.min_config = CircosConfig(self.fasta_dataset, min_rep = 2)
#         self.max_config = CircosConfig(self.fasta_dataset, max_rep = 1)

#         self.min_circos = Circos(self.combined_data, self.min_config)
#         self.max_circos = Circos(self.combined_data, self.max_config)

#         min_path = os.path.join(self.CWD, "circos_test_min.svg")
#         max_path = os.path.join(self.CWD, "circos_test_max.svg")

#         self.min_circos.save(min_path)
#         self.max_circos.save(max_path)

#         min_content = self._read_file(min_path, True)
#         max_content = self._read_file(max_path, True)

#         ref_min_path = os.path.join(self.CWD, "circos_reference_min.svg")
#         ref_max_path = os.path.join(self.CWD, "circos_reference_max.svg")

#         ref_min_content = self._read_file(ref_min_path)
#         ref_max_content = self._read_file(ref_max_path)

#         assert (
#             len(min_content) == len(ref_min_content) and
#             len(max_content) == len(ref_max_content)
#             )

#     def test_positive_save_plot_all_proteins(self):   
#         self.basic_config = CircosConfig(self.fasta_dataset, min_rep = 2)
#         self.modified_config = CircosConfig(self.fasta_dataset, min_rep = 2, plot_all_proteins = True)

#         self.basic_circos = Circos(self.combined_data, self.basic_config)
#         self.modified_circos = Circos(self.combined_data, self.modified_config)

#         basic_path = os.path.join(self.CWD, "circos_test_plot_all_proteins_basic.svg")
#         modified_path = os.path.join(self.CWD, "circos_test_all_proteins_modified.svg")

#         self.basic_circos.save(basic_path)
#         self.modified_circos.save(modified_path)

#         basic_content = self._read_file(basic_path, True)
#         modified_content = self._read_file(modified_path, True)

#         ref_basic_path = os.path.join(self.CWD, "circos_reference_plot_all_proteins_basic.svg")
#         ref_modified_path = os.path.join(self.CWD, "circos_reference_plot_all_proteins_modified.svg")

#         ref_basic_content = self._read_file(ref_basic_path)
#         ref_modified_content = self._read_file(ref_modified_path)

#         assert (
#             len(basic_content) == len(ref_basic_content) and
#             len(modified_content) == len(ref_modified_content)
#             )

#     def test_negative_save_exception_min_and_max_replica(self):   
#         self.config = CircosConfig(self.fasta_dataset, 
#                                    min_rep = 2, 
#                                    max_rep = 1)  # min > max

#         with pytest.raises(ValueError, match = f'ERROR! max_rep is smaller than min_rep'):
#             self.circos = Circos(self.combined_data, self.config)

#     def test_default_set_xls_colors(self):
#         default_intra_color = copy.deepcopy(self.circos.heterotypic_intraprotein_xl_color)
#         default_inter_color = copy.deepcopy(self.circos.heterotypic_interprotein_xl_color)
#         default_homo_color = copy.deepcopy(self.circos.homotypic_xl_color)
#         default_general_color = copy.deepcopy(self.circos.general_xl_color)

#         intra_color = "#21a2ed"
#         inter_color = "#00008B"
#         homo_color = "#ed2b21"
#         general_color = "#7d8082"

#         self.circos.set_xls_colors()

#         assert (
#             # Check default colors
#             default_intra_color == intra_color and
#             default_inter_color == inter_color and
#             default_homo_color == homo_color and
#             default_general_color == general_color and
#             # Check color assignment
#             self.circos.heterotypic_intraprotein_xl_color == intra_color and
#             self.circos.heterotypic_interprotein_xl_color == inter_color and
#             self.circos.homotypic_xl_color == homo_color and
#             self.circos.general_xl_color == general_color
#         )

#     def test_poditive_set_xls_colors(self):
#         light_gray = "#D3D3D3"
#         medium_gray = "#808080"
#         dark_gray = "#404040"
#         charcoal_gray = "#333333"
#         self.circos.set_xls_colors(light_gray, medium_gray, dark_gray, charcoal_gray)

#         assert (
#             self.circos.heterotypic_intraprotein_xl_color == light_gray and
#             self.circos.heterotypic_interprotein_xl_color == medium_gray and
#             self.circos.homotypic_xl_color == dark_gray and
#             self.circos.general_xl_color == charcoal_gray
#         )

#     def test_exception_set_xls_colors(self):
#         light_gray = "Light Grey"
#         medium_gray = "#808080"
#         dark_gray = "#404040"
#         charcoal_gray = "#333333"

#         with pytest.raises(ValueError, match = f'ERROR! Invalid hex color: Light Grey'):
#             self.circos.set_xls_colors(light_gray, medium_gray, dark_gray, charcoal_gray)
