import pytest
import os
import re
from xldg.data import Path, Fasta, MeroX, Domain, CrossLink, ProteinStructure, Util

class TestPath:
    @pytest.fixture(autouse=True)
    def setup(self):
        # Current Working Directory
        self.CWD = os.path.join(os.getcwd(), 'tests', 'files', 'data', 'path')
        self.fasta_files = Path.list_given_type_files(self.CWD, '.fasta')
        self.txt_reference  = [os.path.join(self.CWD, 'abcd1234_file1.txt'),
                               os.path.join(self.CWD, 'abcd1234_file2.txt'),
                               os.path.join(self.CWD, 'abcd1234_file3.txt')]

    def test_positive_list_given_type_files(self):
        sample =  Path.list_given_type_files(self.CWD, '.txt')
        assert sample == self.txt_reference

    def test_negative_list_given_type_files(self):
        sample =  Path.list_given_type_files(self.CWD, '.none')
        reference  = []
        assert sample == reference

    def test_positive_sort_filenames_by_first_integer(self):
        sample = Path.sort_filenames_by_first_integer(self.fasta_files);
        reference  = [os.path.join(self.CWD, 'BSA_1.fasta'),
                      os.path.join(self.CWD, '2_BSA.fasta'),
                      os.path.join(self.CWD, 'BSA_3.fasta')]
        assert sample == reference

    def test_argument_sort_filenames_by_first_integer(self):
        folder_files = Path.list_given_type_files(self.CWD, '.txt')
        original_sample = folder_files.copy()
        folder_files = sorted(folder_files, reverse=True)
        sorted_sample = Path.sort_filenames_by_first_integer(folder_files, 'abcd1234_')

        assert (
        folder_files != original_sample and
        sorted_sample == self.txt_reference
        )

    def test_positive_validate_file_existence(self):
        path = os.path.join(self.CWD, 'BSA_1.fasta')
        sample = Path.validate_file_existence(path)
        # No FileNotFoundError has been raised
        assert 1 == 1

    def test_negative_validate_file_existence(self):
        path = os.path.join(self.CWD, 'None')
        with pytest.raises(FileNotFoundError, match=re.escape(f'File not found: {path}')):
             sample = Path.validate_file_existence(path)

    def test_positive_validate_folder_existence(self):
        sample = Path.validate_folder_existence(self.CWD)
        # No FileNotFoundError has been raised
        assert 1 == 1

    def test_negative_validate_folder_existence(self):
        path = os.path.join(self.CWD, 'None')
        with pytest.raises(FileNotFoundError, match=re.escape(f'Folder not found: {path}')):
             sample = Path.validate_folder_existence(path)

    def test_positive_confirm_file_format(self):
        path = os.path.join(self.CWD, 'BSA_1.fasta')
        fasta_format = 'fasta'
        single_format = Path.confirm_file_format(path, fasta_format)
        many_formats = Path.confirm_file_format(path, 'FASTA', 'txt', 'pdb')

        assert (
            single_format == fasta_format and
            many_formats == fasta_format
            )

    def test_negative_confirm_file_format(self):
        path = os.path.join(self.CWD, 'BSA_1.fasta')
        none_format = 'none'
        with pytest.raises(ValueError, 
                           match=f"Invalid file format. Expected {none_format.upper()}, got FASTA"):
             sample = Path.confirm_file_format(path, none_format)

    def test_read_to_string(self):
        path = os.path.join(self.CWD, 'BSA_1.fasta')
        sample = Path.read_to_string(path)
        reference = ">sp|P02769|ALBU_BOVIN Albumin OS=Bos taurus OX=9913 GN=ALB PE=1 SV=4\n" \
                   "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPF\n" \
                   "DEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEP\n" \
                   "ERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYY\n" \
                   "ANKYNGVFQECCQAEDKGACLLPKIETMREKVLASSARQRLRCASIQKFGERALKAWSVA\n" \
                   "RLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKE\n" \
                   "CCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRR\n" \
                   "HPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEK\n" \
                   "LGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLIL\n" \
                   "NRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLP\n" \
                   "DTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKEACFAVEGPKLVV\n" \
                   "STQTALA"
        assert sample == reference


class TestFasta:
    @pytest.fixture(autouse=True)
    def setup(self):
        # Current Working Directory
        CWD = os.path.join(os.getcwd(), 'tests', 'files', 'data', 'fasta')
        self.single_fasta =  os.path.join(CWD, 'fasta_1.fasta')
        self.all_fasta =  Path.list_given_type_files(CWD, 'fasta')

    def test_single_load_data(self):
        fasta = Fasta.load_data(self.single_fasta, 'Custom')
        assert len(fasta) == 1

    def test_multiple_load_data(self):
        fasta = Fasta.load_data(self.all_fasta, 'Custom')
        assert len(fasta) == 3

class TestDomain:
    @pytest.fixture(autouse=True)
    def setup(self):
        # Current Working Directory
        CWD = os.path.join(os.getcwd(), 'tests', 'files', 'data', 'dmn')
        self.single_domain =  os.path.join(CWD, 'domains_1.dmn')
        self.all_domains =  Path.list_given_type_files(CWD, 'dmn')

    def test_single_load_data(self):
        fasta =Domain.load_data(self.single_domain)
        assert len(fasta) == 10

    def test_multiple_load_data(self):
        fasta = Domain.load_data(self.all_domains)
        assert len(fasta) == 19

class TestMeroX:
    @pytest.fixture(autouse=True)
    def setup(self):
        # Current Working Directory
        CWD = os.path.join(os.getcwd(), 'tests', 'files', 'data', 'merox')
        self.single_zhrm =  os.path.join(CWD, 'dataset_1.zhrm')
        self.all_zhrm =  Path.list_given_type_files(CWD, 'zhrm')

    def test_single_load_data(self):
        crosslinks = MeroX.load_data(self.single_zhrm, 'DSBU')
        assert len(crosslinks) == 227


    def test_multiple_load_data(self):
        crosslinks = MeroX.load_data(self.all_zhrm, 'DSBU')
        assert (
            len(crosslinks) == 3 and
            len(crosslinks[0]) == 227 and
            len(crosslinks[1]) == 88 and
            len(crosslinks[2]) == 112
            )


class TestCrossLink:
    @pytest.fixture(autouse=True)
    def setup(self):
        # Current Working Directory
        self.CWD = os.path.join(os.getcwd(), "tests", "files", "data", "crosslink")
        ZHRM = os.path.join(os.getcwd(), "tests", "files", "data", "merox")
        self.chimerax_folder = os.path.join(self.CWD, 'chimerax')

        zhrm_folder_path = Path.list_given_type_files(ZHRM, 'zhrm')
        self.folder_content = MeroX.load_data(zhrm_folder_path, 'DSBU')
        self.combined_dataset = CrossLink.combine_all(self.folder_content)

    def _read_file(self, file_path: str, delete: bool = False):
        content = None
        with open(file_path, "r", encoding="utf-8") as f:
            content = f.read()

        if delete: 
            os.remove(file_path)
        return content

    def test_negative_filter_by_score(self):
        len_before = len(self.combined_dataset)
        filtered_data = CrossLink.filter_by_score(self.combined_dataset, 0)
        len_after = len(filtered_data)
        assert len_before == len_after

    #TODO Add list argument test
    def test_positive_min_argument_filter_by_score(self):
        filtered_data = CrossLink.filter_by_score(self.combined_dataset, 150)
        assert len(filtered_data) == 25

    def test_positive_max_argument_filter_by_score(self):
        filtered_data = CrossLink.filter_by_score(self.combined_dataset, max_score=1)
        assert len(filtered_data) == 66

    def test_positive_min_and_max_arguments_filter_by_score(self):
        filtered_data = CrossLink.filter_by_score(self.combined_dataset, 120, 150)
        assert len(filtered_data) == 36

    def test_exception_arguments_filter_by_score(self):
        with pytest.raises(ValueError, match = "max_score is smaller than min_score"):
            filtered_data = CrossLink.filter_by_score(self.combined_dataset, 1, 0)

    def test_exception_filter_by_replica(self):
        with pytest.raises(ValueError, match = "max_replica is smaller than min_replica"):
            filtered_data = CrossLink.filter_by_replica(self.combined_dataset, 1, 0)

    def test_negative_filter_by_replica(self):
        len_before = len(self.combined_dataset)
        filtered_data = CrossLink.filter_by_replica(self.combined_dataset)
        len_after = len(filtered_data)
        assert len_before == len_after

    def test_positive_min_argument_filter_by_replica(self):
        filtered_data = CrossLink.filter_by_replica(self.combined_dataset, min_replica=2)
        assert len(filtered_data) == 146

    def test_positive_max_argument_filter_by_replica(self):
        filtered_data = CrossLink.filter_by_replica(self.combined_dataset, max_replica=1)
        assert len(filtered_data) == 281

    def test_positive_min_max_argument_filter_by_replica(self):
        filtered_data = CrossLink.filter_by_replica(self.combined_dataset, 2, 2)
        assert len(filtered_data) == 146
    #TODO Add list argument test to method above

    #TODO Add list argument test to method below
    def test_positive_remove_interprotein_crosslinks(self):
        filtered_data = CrossLink.remove_interprotein(self.combined_dataset)
        assert len(filtered_data) == 296

    def test_positive_remove_intraprotein_crosslinks(self):
        filtered_data = CrossLink.remove_intraprotein(self.combined_dataset)
        assert len(filtered_data) == 138

    def test_positive_remove_homotypic_crosslinks(self):
        filtered_data = CrossLink.remove_homotypic(self.combined_dataset)
        assert len(filtered_data) == 420

    def test_positive_remove_all_xls(self):
        filtered_data = CrossLink.remove_interprotein(self.combined_dataset)
        filtered_data = CrossLink.remove_intraprotein(filtered_data)
        filtered_data = CrossLink.remove_homotypic(filtered_data)
        assert len(filtered_data) == 0

    def test_positive_blank_replica_counter(self):
        len_unmodified_dataset = len(self.combined_dataset)
        filtered_data = CrossLink.filter_by_replica(self.combined_dataset, max_replica=1)
        len_filtered = len(filtered_data)

        len_before = len(self.combined_dataset)
        blanked_data = CrossLink.blank_replica(self.combined_dataset)
        filtered_blank = CrossLink.filter_by_replica(blanked_data, max_replica=1)
        len_after = len(filtered_blank)
        assert (
            len_unmodified_dataset == 427 and 
            len_filtered == 281 and 
            len_before == 427 and
            len_before == len_after
        )

    def test_positive_get_unique(self):
        first_dataset = self.folder_content[0]
        last_dataset = self.folder_content[-1]

        unique_first, unique_last = CrossLink.get_unique(first_dataset, last_dataset)
        common_first, common_last = CrossLink.get_common(first_dataset, last_dataset)
        assert (
            len(unique_first) + len(common_first) == len(first_dataset) and
            len(unique_last) + len(common_last) == len(last_dataset)
        )

    def test_negative_get_unique(self):
        len_before = len(self.combined_dataset)

        reference_dataset = CrossLink.remove_interprotein(self.combined_dataset)
        reference_dataset = CrossLink.remove_intraprotein(reference_dataset)
        reference_dataset = CrossLink.remove_homotypic(reference_dataset)
    
        unique_combined, unique_reference = CrossLink.get_unique(self.combined_dataset, reference_dataset)
        len_after = len(self.combined_dataset)
        assert (
            len_before == len_after and
            len(unique_reference) == 0 and
            len(unique_combined) == len_before
        )

    def test_positive_get_common(self):
        first_dataset = self.folder_content[0] 
        last_dataset = self.folder_content[-1] 

        common_first, common_last = CrossLink.get_common(first_dataset, last_dataset)

        for xl in common_first:
            if xl not in common_last:
                raise ValueError("Common elements are not the same")

        assert len(common_first) == len(common_last)

    def test_negative_get_common(self):
        first_dataset = self.folder_content[0]
        last_dataset = self.folder_content[-1]

        last_dataset = CrossLink.remove_interprotein(last_dataset)
        last_dataset = CrossLink.remove_intraprotein(last_dataset)
        last_dataset = CrossLink.remove_homotypic(last_dataset)

        common_first, common_last = CrossLink.get_common(first_dataset, last_dataset)
        assert len(common_last) == 0 and len(common_first) == 0

    def test_positive_combine_selected(self):
        first_dataset = self.folder_content[0]
        last_dataset = self.folder_content[-1]

        combined_dataset = CrossLink.combine_selected(self.folder_content, [0, 2])
        for xl in combined_dataset:
            if xl in first_dataset or xl in last_dataset:
                continue
            else:
                raise ValueError("Combined elements are not the same")

        assert len(combined_dataset) == 339 and len(combined_dataset) == len(first_dataset) + len(last_dataset)

    def test_negative_combine_selected(self):
        first_dataset = self.folder_content[0]

        last_dataset = CrossLink.remove_interprotein(self.folder_content[-1])
        last_dataset = CrossLink.remove_intraprotein(last_dataset)
        last_dataset = CrossLink.remove_homotypic(last_dataset)
        self.folder_content[2] = last_dataset

        combined_dataset = CrossLink.combine_selected(self.folder_content, [0, 2])
        for xl in combined_dataset:
            if xl not in first_dataset:
                raise ValueError("Combined elements are not the same")

        assert len(combined_dataset) == len(first_dataset)


class TestProteinStructure:
    @pytest.fixture(autouse=True)
    def setup(self):
        # Current Working Directory
        self.CWD = os.path.join(os.getcwd(), 'tests', 'files', 'data', 'structure')

    def test_positive_load_data(self):
        bsa_cif_path = os.path.join(self.CWD, 'BSA.cif')
        bsa_cif_structure = ProteinStructure.load_data(bsa_cif_path)

        bsa_pdb_path = os.path.join(self.CWD, 'BSA.pdb')
        bsa_pdb_structure = ProteinStructure.load_data(bsa_pdb_path)
        assert len(bsa_pdb_structure) == len(bsa_cif_structure) 


class TestUtil:
    def test_positive_generate_list_of_integers(self):
        single_sample = Util.generate_list_of_integers([0, 5])
        multiple_sample = Util.generate_list_of_integers([0, 2], [4, 5])
        
        single_reference = [0, 1, 2, 3, 4, 5]
        multiple_refrence = [0, 1, 2, 4, 5]

        assert (
           single_sample == single_reference and
           multiple_sample == multiple_refrence
        )

    def test_exception_generate_list_of_integers(self):
        with pytest.raises(ValueError, match=f"Start value 5 is greater than end value 0"):
             sample = Util.generate_list_of_integers([5, 0])