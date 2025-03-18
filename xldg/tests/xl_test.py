import pytest
import os
from xldg.utils import PathUtil, DatasetUtil
from xldg.xl import ProteinChainDataset


class TestProteinChainDataset:
    @pytest.fixture(autouse=True)
    def setup(self):
        # Current Working Directory
        self.CWD = os.path.join(os.getcwd(), "tests", "test_data", "xl_test")

    def test_pcd_file_path_argument_constructor_ProteinChainDataset(self):
        path = os.path.join(self.CWD, "non_existent_file.pcd")
        with pytest.raises(FileNotFoundError):
            pcd = ProteinChainDataset(path)

    def test_file_exception_constructor_ProteinChainDataset(self):
        path = os.path.join(self.CWD, "invalid_file.pcd")
        with pytest.raises(ValueError):
            pcd = ProteinChainDataset(path)

    def test_positive_constructor_ProteinChainDataset(self):
        monomer_path = os.path.join(self.CWD, "monomer.pcd")
        monomer_pcd = ProteinChainDataset(monomer_path)
        dimer_path = os.path.join(self.CWD, "dimer.pcd")
        dimer_pcd = ProteinChainDataset(dimer_path)
        assert len(monomer_pcd) == 2 and len(dimer_pcd) == 2


class TestCrossLinkDataset:
    @pytest.fixture(autouse=True)
    def setup(self):
        # Test Data Folder
        TDF = os.path.join(os.getcwd(), "tests", "test_data", "zhrm")
        zhrm_folder_path = PathUtil.list_specified_type_files_from_folder(TDF, '.zhrm')
        folder_content = DatasetUtil.read_merox_zhrm_files_from_path_list(zhrm_folder_path, 'DSBU')
        self.combined_dataset = DatasetUtil.combine_all_datasets(folder_content)

    def test_negative_filter_by_score(self):
        len_before = len(self.combined_dataset)
        self.combined_dataset.filter_by_score(0)
        len_after = len(self.combined_dataset)
        assert len_before == len_after

    def test_positive_min_argument_filter_by_score(self):
        self.combined_dataset.filter_by_score(150)
        assert len(self.combined_dataset) == 25

    def test_positive_max_argument_filter_by_score(self):
        self.combined_dataset.filter_by_score(max_score=1)
        assert len(self.combined_dataset) == 66

    def test_positive_min_and_max_arguments_filter_by_score(self):
        self.combined_dataset.filter_by_score(120, 150)
        assert len(self.combined_dataset) == 36

    def test_exception_arguments_filter_by_score(self):
        with pytest.raises(ValueError, match = "ERROR! max_score is smaller than min_score"):
            filtered_data = DatasetUtil.filter_all_by_score(self.combined_dataset, 1, 0)

    def test_exception_filter_by_replica(self):
        with pytest.raises(ValueError, match = "ERROR! max_rep is smaller than min_rep"):
            self.combined_dataset.filter_by_replica(1, 0)

    def test_negative_filter_by_replica(self):
        len_before = len(self.combined_dataset)
        self.combined_dataset.filter_by_replica()
        len_after = len(self.combined_dataset)
        assert len_before == len_after

    def test_positive_min_argument_filter_by_replica(self):
        self.combined_dataset.filter_by_replica(min_rep=3)
        assert len(self.combined_dataset) == 123

    def test_positive_max_argument_filter_by_replica(self):
        self.combined_dataset.filter_by_replica(max_rep=1)
        assert len(self.combined_dataset) == 246

    def test_positive_min_max_argument_filter_by_replica(self):
        self.combined_dataset.filter_by_replica(2, 2)
        assert len(self.combined_dataset) == 58

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
