import pytest
import os
import copy
from xldg.utils import Path, DatasetUtil
from xldg.xl import ProteinChainDataset, CrossLinkDataset


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
        # Current Working Directory
        self.CWD = os.path.join(os.getcwd(), "tests", "test_data", "xl_test")
        # Test Data Folder
        TDF = os.path.join(os.getcwd(), "tests", "test_data", "zhrm")
        self.chimerax_folder = os.path.join(self.CWD, 'chimerax')

        zhrm_folder_path = Path.list_specified_type_files_from_folder(TDF, '.zhrm')
        self.folder_content = DatasetUtil.read_all_merox_files(zhrm_folder_path, 'DSBU')
        self.combined_dataset = DatasetUtil.combine_all_datasets(self.folder_content)

    def _read_file(self, file_path: str, delete: bool = False):
        content = None
        with open(file_path, "r", encoding="utf-8") as f:
            content = f.read()

        if delete: 
            os.remove(file_path)
        return content

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
        self.combined_dataset.filter_by_replica(min_rep=2)
        assert len(self.combined_dataset) == 146

    def test_positive_max_argument_filter_by_replica(self):
        self.combined_dataset.filter_by_replica(max_rep=1)
        assert len(self.combined_dataset) == 281

    def test_positive_min_max_argument_filter_by_replica(self):
        self.combined_dataset.filter_by_replica(2, 2)
        assert len(self.combined_dataset) == 146

    def test_positive_remove_interprotein_xls(self):
        self.combined_dataset.remove_interprotein_xls()
        assert len(self.combined_dataset) == 296

    def test_positive_remove_intraprotein_xls(self):
        self.combined_dataset.remove_intraprotein_xls()
        assert len(self.combined_dataset) == 138

    def test_positive_remove_homotypic_xls(self):
        self.combined_dataset.remove_homotypic_xls()
        assert len(self.combined_dataset) == 420

    def test_positive_remove_all_xls(self):
        self.combined_dataset.remove_interprotein_xls()
        self.combined_dataset.remove_intraprotein_xls()
        self.combined_dataset.remove_homotypic_xls()
        assert len(self.combined_dataset) == 0

    def test_positive_set_xls_counter_to_one(self):
        unmodified_dataset = self.combined_dataset
        unmodified_dataset.filter_by_replica(max_rep=1)
        len_unmodified_dataset = len(unmodified_dataset)

        len_before = len(self.combined_dataset)
        self.combined_dataset.set_xls_counter_to_one()
        self.combined_dataset.filter_by_replica(max_rep=1)
        len_after = len(self.combined_dataset)
        assert len_unmodified_dataset == 281 and len_before == len_after

    def test_positive_tsv_export_xls_counters(self):
        file_name = 'export_xls_counters.tsv'
        self.combined_dataset.export_xls_counters(self.CWD, file_name)

        content1 = self._read_file(os.path.join(self.CWD, file_name), True)
        content2 = self._read_file(os.path.join(self.CWD, 'counter_reference.tsv'))
        assert content1 == content2

    def test_positive_csv_export_xls_counters(self):
        file_name = 'export_xls_counters.csv'
        self.combined_dataset.export_xls_counters(self.CWD, file_name, ',')

        content1 = self._read_file(os.path.join(self.CWD, file_name), True)
        content2 = self._read_file(os.path.join(self.CWD, 'counter_reference.csv'))
        assert content1 == content2

    def test_positive_monomer_export_for_chimerax(self):
        pcd = ProteinChainDataset(os.path.join(self.CWD, "monomer.pcd"))
        self.combined_dataset.set_xls_counter_to_one()
        self.combined_dataset.remove_interprotein_xls()
        self.combined_dataset.export_for_chimerax(pcd, self.chimerax_folder, "monomer")

        content1 = self._read_file(os.path.join(self.chimerax_folder, "monomer_heterotypical_intraprotein_xl_1_rep.pb"), True)
        content2 = self._read_file(os.path.join(self.chimerax_folder, 'monomer_reference.pb'))
        assert content1 == content2

    def test_positive_dimer_export_for_chimerax(self):
        pcd = ProteinChainDataset(os.path.join(self.CWD, "dimer.pcd"))
        self.combined_dataset.set_xls_counter_to_one()
        self.combined_dataset.remove_intraprotein_xls()
        self.combined_dataset.remove_homotypic_xls()
        self.combined_dataset.export_for_chimerax(pcd, self.chimerax_folder, "dimer")

        content1 = self._read_file(os.path.join(self.chimerax_folder, "dimer_heterotypical_interprotein_xl_1_rep.pb"), True)
        content2 = self._read_file(os.path.join(self.chimerax_folder, 'dimer_reference.pb'))
        assert content1 == content2

    def test_positive_arguments_export_for_chimerax(self):
        pcd = ProteinChainDataset(os.path.join(self.CWD, "dimer.pcd"))
        self.combined_dataset.set_xls_counter_to_one()
        self.combined_dataset.export_for_chimerax(pcd, self.chimerax_folder, "color", 0.3, 10, "#D3D3D3", "#808080", "#404040")

        interprotein_pb = self._read_file(os.path.join(self.chimerax_folder, "color_heterotypical_interprotein_xl_1_rep.pb"), True)
        intraprotein_pb = self._read_file(os.path.join(self.chimerax_folder, "color_heterotypical_intraprotein_xl_1_rep.pb"), True)
        homotypic_pb = self._read_file(os.path.join(self.chimerax_folder, "color_homotypical_xl_1_rep.pb"), True)

        ref_interprotein_pb = self._read_file(os.path.join(self.chimerax_folder, "color_heterotypical_interprotein_reference.pb"))
        ref_intraprotein_pb = self._read_file(os.path.join(self.chimerax_folder, "color_heterotypical_intraprotein_reference.pb"))
        ref_homotypic_pb = self._read_file(os.path.join(self.chimerax_folder, "color_homotypical_reference.pb"))
        assert (
            interprotein_pb == ref_interprotein_pb and
            intraprotein_pb == ref_intraprotein_pb and
            homotypic_pb == ref_homotypic_pb
        )

    def test_positive_export_ppis_for_gephi(self):
        pcd_monomer = ProteinChainDataset(os.path.join(self.CWD, "monomer.pcd"))
        pcd_dimer = ProteinChainDataset(os.path.join(self.CWD, "dimer.pcd"))
        
        save_monomer = "ppis_for_gephi_test_monomer.gexf"
        save_dimer = "ppis_for_gephi_test_dimer.gexf"

        self.combined_dataset.export_ppis_for_gephi(pcd_monomer, self.chimerax_folder, save_monomer)
        self.combined_dataset.export_ppis_for_gephi(pcd_dimer, self.chimerax_folder, save_dimer)

        monomer_sample_path = os.path.join(self.chimerax_folder, save_monomer)
        monomer_reference_path = os.path.join(self.chimerax_folder, "ppis_for_gephi_reference_monomer.gexf")

        dimer_sample_path = os.path.join(self.chimerax_folder, save_dimer)
        dimer_reference_path = os.path.join(self.chimerax_folder, "ppis_for_gephi_reference_dimer.gexf")

        monomer_sample = self._read_file(monomer_sample_path, True)
        monomer_reference = self._read_file(monomer_reference_path)

        dimer_sample = self._read_file(dimer_sample_path, True)
        dimer_reference = self._read_file(dimer_reference_path)

        assert (
           dimer_sample == dimer_reference and
           monomer_sample == monomer_reference
        )

    def test_exception_export_ppis_for_gephi(self):
        pcd = ProteinChainDataset(os.path.join(self.CWD, "dimer.pcd"))
        save_name = "ppis_for_gephi_test.PDF"
        
        with pytest.raises(ValueError, match = f'ERROR! Wrong file extension in {save_name}. Only ".gexf" format is supported'):
            self.combined_dataset.export_ppis_for_gephi(pcd, self.chimerax_folder, save_name)

    def test_positive_export_aais_for_gephi(self):
        pcd_monomer = ProteinChainDataset(os.path.join(self.CWD, "monomer.pcd"))
        pcd_dimer = ProteinChainDataset(os.path.join(self.CWD, "dimer.pcd"))
        
        save_monomer = "aais_for_gephi_test_monomer.gexf"
        save_dimer = "aais_for_gephi_test_dimer.gexf"

        self.combined_dataset.export_aais_for_gephi(pcd_monomer, self.chimerax_folder, save_monomer)
        self.combined_dataset.export_aais_for_gephi(pcd_dimer, self.chimerax_folder, save_dimer)

        monomer_sample_path = os.path.join(self.chimerax_folder, save_monomer)
        monomer_reference_path = os.path.join(self.chimerax_folder, "aais_for_gephi_reference_monomer.gexf")

        dimer_sample_path = os.path.join(self.chimerax_folder, save_dimer)
        dimer_reference_path = os.path.join(self.chimerax_folder, "aais_for_gephi_reference_dimer.gexf")

        monomer_sample = self._read_file(monomer_sample_path, True)
        monomer_reference = self._read_file(monomer_reference_path)

        dimer_sample = self._read_file(dimer_sample_path, True)
        dimer_reference = self._read_file(dimer_reference_path)

        assert (
           dimer_sample == dimer_reference and
           monomer_sample == monomer_reference
        )

    def test_exception_export_aais_for_gephi(self):
        pcd = ProteinChainDataset(os.path.join(self.CWD, "dimer.pcd"))
        save_name = "ppis_for_gephi_test.PDF"
        
        with pytest.raises(ValueError, match = f'ERROR! Wrong file extension in {save_name}. Only ".gexf" format is supported'):
            self.combined_dataset.export_aais_for_gephi(pcd, self.chimerax_folder, save_name)

    def test_positive_unique_elements(self):
        first_dataset = self.folder_content[0]
        last_dataset = self.folder_content[-1]

        unique_first, unique_last = CrossLinkDataset.unique_elements(first_dataset, last_dataset)
        common_first, common_last = CrossLinkDataset.common_elements(first_dataset, last_dataset)
        assert (
            len(unique_first) + len(common_first) == len(first_dataset) and
            len(unique_last) + len(common_last) == len(last_dataset)
        )

    def test_negative_unique_elements(self):
        len_before = len(self.combined_dataset)

        reference_dataset = copy.deepcopy(self.combined_dataset)
        reference_dataset.remove_intraprotein_xls()
        reference_dataset.remove_interprotein_xls()
        reference_dataset.remove_homotypic_xls()
    
        unique_combined, unique_reference = CrossLinkDataset.unique_elements(self.combined_dataset, reference_dataset)
        len_after = len(self.combined_dataset)
        assert (
            len_before == len_after and
            len(unique_reference) == 0 and
            len(unique_combined) == len_before
        )

    def test_positive_common_elements(self):
        first_dataset = self.folder_content[0]
        last_dataset = self.folder_content[-1]

        common_first, common_last = CrossLinkDataset.common_elements(first_dataset, last_dataset)

        for xl in common_first:
            if xl not in common_last:
                raise ValueError("Common elements are not the same")

        assert len(common_first) == len(common_last)

    def test_negative_common_elements(self):
        first_dataset = self.folder_content[0]
        last_dataset = self.folder_content[-1]

        last_dataset.remove_intraprotein_xls()
        last_dataset.remove_interprotein_xls()
        last_dataset.remove_homotypic_xls()

        common_first, common_last = CrossLinkDataset.common_elements(first_dataset, last_dataset)
        for xl in common_first:
            if xl not in first_dataset:
                raise ValueError("Common elements are not the same")

        assert len(common_last) == 0 and len(common_first) == 0

    def test_positive_combine_datasets(self):
        first_dataset = self.folder_content[0]
        last_dataset = self.folder_content[-1]

        combined_dataset = CrossLinkDataset.combine_datasets([first_dataset, last_dataset])
        for xl in combined_dataset:
            if xl in first_dataset or xl in last_dataset:
                continue
            else:
                raise ValueError("Combined elements are not the same")

        assert len(combined_dataset) == 339 and len(combined_dataset) == len(first_dataset) + len(last_dataset)

    def test_negative_combine_datasets(self):
        first_dataset = self.folder_content[0]
        last_dataset = self.folder_content[-1]

        last_dataset.remove_intraprotein_xls()
        last_dataset.remove_interprotein_xls()
        last_dataset.remove_homotypic_xls()

        combined_dataset = CrossLinkDataset.combine_datasets([first_dataset, last_dataset])
        for xl in combined_dataset:
            if xl not in first_dataset:
                raise ValueError("Combined elements are not the same")

        assert len(combined_dataset) == len(first_dataset)
