import pytest
import os
import sys

from xldg.utils import PathUtil, DatasetUtil
from xldg.fasta import Fasta, FastaDataset


class TestFasta:
    @pytest.fixture(autouse=True)
    def setup(self):
        # Current Working Directory
        CWD = os.path.join(os.getcwd(), "tests", "test_data", "fasta")

    def test_positive_Uniprot_constructor(self):
        header = '>sp|P02769|ALBU_BOVIN Albumin OS=Bos taurus OX=9913 GN=ALB PE=1 SV=4'
        sequence = 'MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVN'
        fasta_format = 'Uniprot'
        remove_parenthesis = False
        fasta = Fasta(header, sequence, fasta_format, remove_parenthesis)
        assert fasta.db_id == 'P02769' and fasta.prot_gene == 'ALB' and fasta.seq_length == len(sequence) + 2

    def test_positive_Araport11_constructor(self):
        header = '>AT4G35310.1 | Symbols: CPK5, ATCPK5 | calmodulin-domain protein kinase 5 | chr4:16802436-16804628 FORWARD LENGTH=556'
        sequence = 'MGNSCRGSFKDKLDEGDNNKPEDYSKTSTTNLSSNSDHSPNAADIIAQEFSKDNNSNNNSKDPALVIPLREPIMRRNPDN'
        fasta_format = 'Araport11'
        remove_parenthesis = False
        fasta = Fasta(header, sequence, fasta_format, remove_parenthesis)
        assert fasta.db_id == 'AT4G35310.1' and fasta.prot_gene == 'CPK5' and fasta.seq_length == len(sequence) + 2

    def test_positive_Custom_constructor(self):
        header = '>BSA'
        sequence = 'MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVN'
        fasta_format = 'Custom'
        remove_parenthesis = False
        fasta = Fasta(header, sequence, fasta_format, remove_parenthesis)
        assert (
            fasta.db_id == header and 
            fasta.prot_gene == header and 
            fasta.seq_length == len(sequence) + 2
        )

    def test_remove_parenthesis_constructor(self):
        original_header = '>sp|P02769|(ALBU_BOVIN) (Albumin) OS=(Bos taurus) OX=(9913) GN=(ALB) PE=(1) SV=(4)'
        clean_header = '>sp|P02769|ALBU_BOVIN Albumin OS=Bos taurus OX=9913 GN=ALB PE=1 SV=4'
        sequence = 'M'
        fasta_format = 'Custom'
        remove_parenthesis = True
        fasta = Fasta(original_header, sequence, fasta_format, remove_parenthesis)
        assert fasta.db_id == clean_header and fasta.prot_gene == clean_header

    def test_nonexistent_format_exception_constructor(self):
        header = '>AT4G35310.1 | Symbols: CPK5, ATCPK5 | calmodulin-domain protein kinase 5 | chr4:16802436-16804628 FORWARD LENGTH=556'
        sequence = 'M'
        wrong_format = 'Nonexistent'
        remove_parenthesis = False
        with pytest.raises(ValueError, match = "ERROR! Wrong FASTA format: Nonexistent"):
            fasta = Fasta(header, sequence, wrong_format, remove_parenthesis)

    def test_wrong_format_exception_constructor_Fasta(self):
        header = '>CPK5'
        sequence = 'M'
        wrong_format = 'Uniprot'
        remove_parenthesis = False
        with pytest.raises(ValueError, match = "ERROR! Wrong FASTA format: Uniprot"):
            fasta = Fasta(header, sequence, wrong_format, remove_parenthesis)


class TestFastaDataset:
    @pytest.fixture(autouse=True)
    def setup(self):
        # Current Working Directory
        self.CWD = os.path.join(os.getcwd(), "tests", "test_data", "fasta")
        # Test Data Folder
        TDF = os.path.join(os.getcwd(), "tests", "test_data", "zhrm")

        self.fasta_files = PathUtil.list_specified_type_files_from_folder(self.CWD, ".fasta")
        self.fasta_dataset = FastaDataset(self.fasta_files, "Custom", False)

        zhrm_folder_path = PathUtil.list_specified_type_files_from_folder(TDF, '.zhrm')
        folder_content = DatasetUtil.read_merox_zhrm_files_from_path_list(zhrm_folder_path, 'DSBU')
        self.combined_replicas = DatasetUtil.combine_all_datasets(folder_content)

        self.fas_files = PathUtil.list_specified_type_files_from_folder(self.CWD, ".fas")
        self.fas_dataset = FastaDataset(self.fas_files, "Uniprot")

    def test_positive_constructor(self):
        fasta_dataset = FastaDataset(self.fasta_files, "Custom")
        assert len(fasta_dataset) == 3

    def test_nonexistent_format_exception_constructor(self):
        with pytest.raises(ValueError, match = "ERROR! Wrong FASTA format: Nonexistent"):
            fasta_dataset = FastaDataset(self.fasta_files, "Nonexistent")

    def test_wrong_format_exception_constructor(self):
        with pytest.raises(ValueError, match = "ERROR! Wrong FASTA format: Uniprot"):
            # Custom format is expected
            fasta_dataset = FastaDataset(self.fasta_files, "Uniprot")

    def test_positive_filter_by_CrossLinkDataset(self):
        fasta_dataset = FastaDataset(self.fasta_files, "Custom", False)
        fasta_dataset.filter_by_CrossLinkDataset(self.combined_replicas)
        assert len(fasta_dataset) == 2

    def test_negative_filter_by_CrossLinkDataset(self):
        self.combined_replicas.filter_by_score(min_score = sys.maxsize - 1, 
                                          max_score = sys.maxsize)

        self.fasta_dataset.filter_by_CrossLinkDataset(self.combined_replicas)
        assert len(self.fasta_dataset) == 0

    def test_positive_find_gene_by_fasta_header(self):
        header = r'>sp|P02358|RS6_ECOLI Small ribosomal subunit protein bS6 OS=Escherichia coli (strain K12) OX=83333 GN=rpsF PE=1 SV=1'
        fasta_match = self.fas_dataset.find_gene_by_fasta_header(header)
        expected_match = 'rpsF'
        assert fasta_match == expected_match

    def test_negative_find_gene_by_fasta_header(self):
        header = r'None'
        fasta_match = self.fas_dataset.find_gene_by_fasta_header(header)
        assert fasta_match == None

    def test_positive_remove_parenthesis_argument_constructor(self):
        fas_dataset = FastaDataset(self.fas_files, "Uniprot", True)
        header = r'>sp|P02358|RS6_ECOLI Small ribosomal subunit protein bS6 OS=Escherichia coli strain K12 OX=83333 GN=rpsF PE=1 SV=1'
        fasta_match = fas_dataset.find_gene_by_fasta_header(header)
        expected_match = 'rpsF'
        assert fasta_match == expected_match

    def test_save(self):
        self.fasta_dataset.filter_by_CrossLinkDataset(self.combined_replicas)

        saved_fasta = os.path.join(self.CWD, "test.fasta")
        if os.path.exists(saved_fasta):
            raise FileExistsError("Fasta test file should not exist")

        self.fasta_dataset.save(self.CWD, "test.fasta")

        if not os.path.exists(saved_fasta):
            raise FileExistsError("Fasta test file should exists")

        saved_fasta_dataset = FastaDataset([saved_fasta], "Custom", False)
        os.remove(saved_fasta)
        assert str(saved_fasta_dataset) == str(self.fasta_dataset)
