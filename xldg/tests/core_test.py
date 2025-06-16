import pytest
import os
import sys
import copy
from typing import List

from xldg.core import Node, Atom, ProteinStructureDataset, FastaEntity, FastaDataset, DomainEntity, DomainDataset, ProteinChainDataset, CrossLinkDataset
from xldg.data import Path, MeroX, CrossLink


class TestProteinChainDataset:
    @pytest.fixture(autouse=True)
    def setup(self):
        # Current Working Directory
        self.CWD = os.path.join(os.getcwd(), 'tests', 'files', 'data', 'pcd')

    def test_file_exception_constructor(self):
        path = os.path.join(self.CWD, 'invalid_file.pcd')
        with pytest.raises(ValueError):
            pcd = ProteinChainDataset(path)

    def test_positive_constructor_single_chain(self):
        monomer_path = os.path.join(self.CWD, 'monomer.pcd')
        monomer_content = Path.read_to_string(monomer_path)
        monomer_pcd = ProteinChainDataset( monomer_content)
        assert len(monomer_pcd) == 2

    def test_positive_constructor_multiple_chains(self):
        dimer_path = os.path.join(self.CWD, 'dimer.pcd')
        dimer_content = Path.read_to_string(dimer_path)
        dimer_pcd = ProteinChainDataset(dimer_content)
        assert len(dimer_pcd) == 2


class TestNode:
    def test_distance_to(self):
        node1 = Node(1.0, 2.0, 3.0)
        node2 = Node(4.0, 6.0, 8.0)
        distance = node1.distance_to(node2)
        assert  round(distance, 3) == 7.071

    def test_to_tuple(self):
        node = Node(1.0, 2.0, 3.0)
        to_tuple = node.to_tuple()
        assert to_tuple == (1.0, 2.0, 3.0)


class TestAtom:
    def test_distance_to(self):
        atom1 = Atom(1, '12A', 'CA', 'A', 1.0, 2.0, 3.0)
        atom2 = Atom(2, '15B', 'CB', 'A', 4.0, 6.0, 8.0)
        distance = atom1.distance_to(atom2)
        assert  round(distance, 3) == 7.071


class TestProteinStructureDataset:
    @pytest.fixture(autouse=True)
    def setup(self):
        # Current Working Directory
        self.CWD = os.path.join(os.getcwd(), 'tests', 'files', 'data', 'structure')

        pcd_path = os.path.join(os.getcwd(), 'tests', 'files', 'data', 'pcd', 'monomer.pcd')
        content = Path.read_to_string(pcd_path)
        self.pcd = ProteinChainDataset(content)

        self.structure_path = os.path.join(self.CWD, 'two_monomers_complex_structure_reference.cif')
        structure_content = Path.read_to_string(self.structure_path)
        self.structure = ProteinStructureDataset(structure_content, 'cif')

    def test_positive_constructor(self):
        cif_path = os.path.join(self.CWD, 'BSA.cif')
        cif_content = Path.read_to_string(cif_path)
        cif_structure = ProteinStructureDataset(cif_content, 'cif')

        pdb_path = os.path.join(self.CWD, 'BSA.pdb')
        pdb_content = Path.read_to_string(pdb_path)
        pdb_structure = ProteinStructureDataset(pdb_content, 'pdb')
        assert (
            len(pdb_structure) == len(cif_structure) and 
            len(pdb_structure) == 4854
            )

    def test_negative_constructor(self):
        bsa_cif_path = os.path.join(self.CWD, 'BSA.cif')
        with pytest.raises(ValueError, match = 'Unsupported file format. Only .pdb and .cif formats can be used'):
            wrong_type = ProteinStructureDataset(bsa_cif_path, 'none')

    def test_exception_wrong_atom_type_argument_predict_crosslinks(self):
        with pytest.raises(ValueError, match = f'Invalid atom type ?'):
            self.structure.predict_crosslinks(self.pcd, '{K', '{K', atom_type='?')

    def test_exception_node_multiplier_argument_predict_crosslinks(self):
        with pytest.raises(ValueError, match = 'Node multiplier must be greater than 0, got 0'):
            self.structure.predict_crosslinks(self.pcd, '{K', '{K', node_multiplier=0)

    def test_exception_min_length_argument_predict_crosslinks(self):
        with pytest.raises(ValueError, match = 'Min length 0 must be greater than 0'):
            self.structure.predict_crosslinks(self.pcd, '{K', '{K', min_length=0)

    def test_exception_min_and_max_length_argument_predict_crosslinks(self):
        with pytest.raises(ValueError, match = 'Max length 5 must be greater than min length 10'):
            self.structure.predict_crosslinks(self.pcd, '{K', '{K', min_length=10, max_length=5)

    def test_positive_direct_path_predict_crosslinks(self):
        prediction = self.structure.predict_crosslinks(self.pcd, 'M', 'M', min_length=10, max_length=35, direct_path=True, num_processes=0)
        assert len(prediction) == 154

    def test_positive_indirect_path_predict_crosslinks(self):
        prediction = self.structure.predict_crosslinks(self.pcd, 'M', 'M', min_length=10, max_length=35, direct_path=False, num_processes=0)
        assert len(prediction) == 126


class TestCrossLinkDataset:
    @pytest.fixture(autouse=True)
    def setup(self):
        # Current Working Directory
        self.CWD = os.path.join(os.getcwd(), 'tests', 'files', 'data', 'crosslink')
        self.pcd = os.path.join(os.getcwd(), 'tests', 'files', 'data', 'pcd')
        # Test Data Folder
        TDF = os.path.join(os.getcwd(), 'tests', 'files', 'data', 'merox')

        self.chimerax_folder = os.path.join(self.CWD, 'chimerax')
        self.gephi_folder = os.path.join(self.CWD, 'gephi')

        zhrm_folder_path = Path.list_given_type_files(TDF, 'zhrm')
        self.folder_content = MeroX.load_data(zhrm_folder_path, 'DSBU')
        self.combined_dataset = CrossLink.combine_all(self.folder_content)
        pcd_content = Path.read_to_string(os.path.join(self.pcd, 'dimer.pcd'))
        self.dimer_pcd = ProteinChainDataset(pcd_content)

    def _read_file(self, file_path: str, delete: bool = False):
        content = None
        with open(file_path, 'r', encoding='utf-8') as f:
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

    def test_positive_remove_interprotein_crosslinks(self):
        self.combined_dataset.remove_interprotein_crosslinks()
        assert len(self.combined_dataset) == 296

    def test_positive_remove_intraprotein_crosslinks(self):
        self.combined_dataset.remove_intraprotein_crosslinks()
        assert len(self.combined_dataset) == 138

    def test_positive_remove_homeotypic_crosslinks(self):
        self.combined_dataset.remove_homeotypic_crosslinks()
        assert len(self.combined_dataset) == 420

    def test_positive_remove_all_xls(self):
        self.combined_dataset.remove_interprotein_crosslinks()
        self.combined_dataset.remove_intraprotein_crosslinks()
        self.combined_dataset.remove_homeotypic_crosslinks()
        assert len(self.combined_dataset) == 0

    def test_positive_blank_replica_counter(self):
        unmodified_dataset = self.combined_dataset
        unmodified_dataset.filter_by_replica(max_rep=1)
        len_unmodified_dataset = len(unmodified_dataset)

        len_before = len(self.combined_dataset)
        self.combined_dataset.blank_replica_counter()
        self.combined_dataset.filter_by_replica(max_rep=1)
        len_after = len(self.combined_dataset)
        assert len_unmodified_dataset == 281 and len_before == len_after

    def test_positive_tsv_save_crosslink_counter(self):
        file_name = 'save_crosslink_counter.tsv'
        self.combined_dataset.save_crosslink_counter(self.CWD, file_name)

        content1 = self._read_file(os.path.join(self.CWD, file_name), True)
        content2 = self._read_file(os.path.join(self.CWD, 'counter_reference.tsv'))
        assert content1 == content2

    def test_positive_csv_save_crosslink_counter(self):
        file_name = 'save_crosslink_counter.csv'
        self.combined_dataset.save_crosslink_counter(self.CWD, file_name, ',')

        content1 = self._read_file(os.path.join(self.CWD, file_name), True)
        content2 = self._read_file(os.path.join(self.CWD, 'counter_reference.csv'))
        assert content1 == content2

    def test_positive_monomer_export_for_chimerax(self):
        pcd_content = Path.read_to_string(os.path.join(self.pcd, 'monomer.pcd'))
        pcd = ProteinChainDataset(pcd_content)
        self.combined_dataset.blank_replica_counter()
        self.combined_dataset.remove_interprotein_crosslinks()
        self.combined_dataset.export_for_chimerax(pcd, self.chimerax_folder, 'monomer')

        content1 = self._read_file(os.path.join(self.chimerax_folder, 'monomer_intraprotein_xl_1_rep.pb'), True)
        content2 = self._read_file(os.path.join(self.chimerax_folder, 'monomer_reference.pb'))
        assert len(content1) == len(content2)

    def test_positive_dimer_export_for_chimerax(self):
        self.combined_dataset.blank_replica_counter()
        self.combined_dataset.remove_intraprotein_crosslinks()
        self.combined_dataset.remove_homeotypic_crosslinks()
        self.combined_dataset.export_for_chimerax(self.dimer_pcd, self.chimerax_folder, 'dimer')

        content1 = self._read_file(os.path.join(self.chimerax_folder, 'dimer_interprotein_xl_1_rep.pb'), True)
        content2 = self._read_file(os.path.join(self.chimerax_folder, 'dimer_reference.pb'))
        assert len(content1) == len(content2)

    def test_positive_cosmetics_export_for_chimerax(self):
        self.combined_dataset.blank_replica_counter()
        self.combined_dataset.export_for_chimerax(self.dimer_pcd, self.chimerax_folder, 'color', 0.3, 10, '#D3D3D3', '#808080')

        interprotein_pb = self._read_file(os.path.join(self.chimerax_folder, 'color_interprotein_xl_1_rep.pb'), True)
        intraprotein_pb = self._read_file(os.path.join(self.chimerax_folder, 'color_intraprotein_xl_1_rep.pb'), True)
        homeotypic_pb = self._read_file(os.path.join(self.chimerax_folder, 'color_homeotypical_xl_1_rep.pb'), True)

        ref_interprotein_pb = self._read_file(os.path.join(self.chimerax_folder, 'color_interprotein_reference.pb'))
        ref_intraprotein_pb = self._read_file(os.path.join(self.chimerax_folder, 'color_intraprotein_reference.pb'))
        ref_homeotypic_pb = self._read_file(os.path.join(self.chimerax_folder, 'color_homeotypical_reference.pb'))
        assert(
            len(interprotein_pb) == len(ref_interprotein_pb) and
            len(intraprotein_pb) == len(ref_intraprotein_pb) and
            len(homeotypic_pb) == len(ref_homeotypic_pb)
        )

    def test_positive_export_ppis_for_gephi(self):
        monomer_content = Path.read_to_string(os.path.join(self.pcd, 'monomer.pcd'))
        dimer_content = Path.read_to_string(os.path.join(self.pcd, 'dimer.pcd'))

        pcd_monomer = ProteinChainDataset(monomer_content)
        pcd_dimer = ProteinChainDataset(dimer_content)
        
        save_monomer = 'ppis_for_gephi_test_monomer.gexf'
        save_dimer = 'ppis_for_gephi_test_dimer.gexf'

        self.combined_dataset.export_ppis_for_gephi(self.gephi_folder, save_monomer, pcd_monomer)
        self.combined_dataset.export_ppis_for_gephi(self.gephi_folder, save_dimer, pcd_dimer)

        monomer_sample_path = os.path.join(self.gephi_folder, save_monomer)
        monomer_reference_path = os.path.join(self.gephi_folder, 'ppis_for_gephi_reference_monomer.gexf')

        dimer_sample_path = os.path.join(self.gephi_folder, save_dimer)
        dimer_reference_path = os.path.join(self.gephi_folder, 'ppis_for_gephi_reference_dimer.gexf')

        monomer_sample = self._read_file(monomer_sample_path, True)
        monomer_reference = self._read_file(monomer_reference_path)

        dimer_sample = self._read_file(dimer_sample_path, True)
        dimer_reference = self._read_file(dimer_reference_path)

        assert(
           len(dimer_sample) == len(dimer_reference) and
           len(monomer_sample) == len(monomer_reference)
        )

    def test_exception_export_ppis_for_gephi(self):
        pcd_content = Path.read_to_string(os.path.join(self.pcd, 'dimer.pcd'))
        pcd = ProteinChainDataset(pcd_content)
        save_name = 'ppis_for_gephi_test.PDF'
        
        with pytest.raises(ValueError, match = f'Wrong data format is provided in {save_name}. Only ".gexf" format is supported'):
            self.combined_dataset.export_ppis_for_gephi(self.chimerax_folder, save_name, pcd)

    def test_positive_export_aais_for_gephi(self):
        monomer_content = Path.read_to_string(os.path.join(self.pcd, 'monomer.pcd'))
        dimer_content = Path.read_to_string(os.path.join(self.pcd, 'dimer.pcd'))
        
        pcd_monomer = ProteinChainDataset(monomer_content)
        pcd_dimer = ProteinChainDataset(dimer_content)
        
        save_monomer = 'aais_for_gephi_test_monomer.gexf'
        save_dimer = 'aais_for_gephi_test_dimer.gexf'
        
        self.combined_dataset.export_aais_for_gephi(pcd_monomer, self.gephi_folder, save_monomer)
        self.combined_dataset.export_aais_for_gephi(pcd_dimer, self.gephi_folder, save_dimer)

        monomer_sample_path = os.path.join(self.gephi_folder, save_monomer)
        monomer_reference_path = os.path.join(self.gephi_folder, 'aais_for_gephi_reference_monomer.gexf')

        dimer_sample_path = os.path.join(self.gephi_folder, save_dimer)
        dimer_reference_path = os.path.join(self.gephi_folder, 'aais_for_gephi_reference_dimer.gexf')

        monomer_sample = self._read_file(monomer_sample_path, True)
        monomer_reference = self._read_file(monomer_reference_path)

        dimer_sample = self._read_file(dimer_sample_path, True)
        dimer_reference = self._read_file(dimer_reference_path)

        assert(
           dimer_sample == dimer_reference and
           monomer_sample == monomer_reference
        )

    def test_exception_export_aais_for_gephi(self):
        pcd_content = Path.read_to_string(os.path.join(self.pcd, 'dimer.pcd'))
        pcd = ProteinChainDataset(pcd_content)
        save_name = 'ppis_for_gephi_test.PDF'
        
        with pytest.raises(ValueError, match = f'Wrong data format is provided in {save_name}. Only ".gexf" format is supported'):
            self.combined_dataset.export_aais_for_gephi(pcd, self.gephi_folder, save_name)

    def test_positive_unique_elements(self):
        first_dataset = self.folder_content[0]
        last_dataset = self.folder_content[-1]

        unique_first, unique_last = CrossLinkDataset.unique_elements(first_dataset, last_dataset)
        common_first, common_last = CrossLinkDataset.common_elements(first_dataset, last_dataset)
        assert(
            len(unique_first) + len(common_first) == len(first_dataset) and
            len(unique_last) + len(common_last) == len(last_dataset)
        )

    def test_negative_unique_elements(self):
        len_before = len(self.combined_dataset)

        reference_dataset = copy.deepcopy(self.combined_dataset)
        reference_dataset.remove_intraprotein_crosslinks()
        reference_dataset.remove_interprotein_crosslinks()
        reference_dataset.remove_homeotypic_crosslinks()
    
        unique_combined, unique_reference = CrossLinkDataset.unique_elements(self.combined_dataset, reference_dataset)
        len_after = len(self.combined_dataset)
        assert(
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
                raise ValueError('Common elements are not the same')

        assert len(common_first) == len(common_last)

    def test_negative_common_elements(self):
        first_dataset = self.folder_content[0]
        last_dataset = self.folder_content[-1]

        last_dataset.remove_intraprotein_crosslinks()
        last_dataset.remove_interprotein_crosslinks()
        last_dataset.remove_homeotypic_crosslinks()

        common_first, common_last = CrossLinkDataset.common_elements(first_dataset, last_dataset)
        for xl in common_first:
            if xl not in first_dataset:
                raise ValueError('Common elements are not the same')

        assert len(common_last) == 0 and len(common_first) == 0

    def test_positive_combine_datasets(self):
        first_dataset = self.folder_content[0]
        last_dataset = self.folder_content[-1]

        combined_dataset = first_dataset + last_dataset
        for xl in combined_dataset:
            if xl in first_dataset or xl in last_dataset:
                continue
            else:
                raise ValueError('Combined elements are not the same')

        assert len(combined_dataset) == 339 and len(combined_dataset) == len(first_dataset) + len(last_dataset)


class TestFastaEntity:
    @pytest.fixture(autouse=True)
    def setup(self):
        # Current Working Directory
        CWD = os.path.join(os.getcwd(), 'tests', 'test_data', 'fasta')

    def test_positive_Uniprot_constructor(self):
        header = '>sp|P02769|ALBU_BOVIN Albumin OS=Bos taurus OX=9913 GN=ALB PE=1 SV=4'
        sequence = 'MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVN'
        fasta_format = 'Uniprot'
        remove_parenthesis = False
        fasta = FastaEntity(header, sequence, fasta_format, remove_parenthesis)
        assert fasta.db_id == 'P02769' and fasta.prot_gene == 'ALB' and fasta.seq_length == len(sequence) + 2

    def test_positive_Araport11_constructor(self):
        header = '>AT4G35310.1 | Symbols: CPK5, ATCPK5 | calmodulin-domain protein kinase 5 | chr4:16802436-16804628 FORWARD LENGTH=556'
        sequence = 'MGNSCRGSFKDKLDEGDNNKPEDYSKTSTTNLSSNSDHSPNAADIIAQEFSKDNNSNNNSKDPALVIPLREPIMRRNPDN'
        fasta_format = 'Araport11'
        remove_parenthesis = False
        fasta = FastaEntity(header, sequence, fasta_format, remove_parenthesis)
        assert fasta.db_id == 'AT4G35310.1' and fasta.prot_gene == 'CPK5' and fasta.seq_length == len(sequence) + 2

    def test_positive_Custom_constructor(self):
        header = '>BSA'
        sequence = 'MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVN'
        fasta_format = 'Custom'
        remove_parenthesis = False
        fasta = FastaEntity(header, sequence, fasta_format, remove_parenthesis)
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
        fasta = FastaEntity(original_header, sequence, fasta_format, remove_parenthesis)
        assert fasta.db_id == clean_header and fasta.prot_gene == clean_header

    def test_nonexistent_format_exception_constructor(self):
        header = '>AT4G35310.1 | Symbols: CPK5, ATCPK5 | calmodulin-domain protein kinase 5 | chr4:16802436-16804628 FORWARD LENGTH=556'
        sequence = 'M'
        wrong_format = 'Nonexistent'
        remove_parenthesis = False
        with pytest.raises(ValueError, match = 'Wrong FASTA format: Nonexistent'):
            fasta = FastaEntity(header, sequence, wrong_format, remove_parenthesis)

    def test_wrong_format_exception_constructor_FastaEntity(self):
        header = '>CPK5'
        sequence = 'M'
        wrong_format = 'Uniprot'
        remove_parenthesis = False
        with pytest.raises(ValueError, match = 'Wrong FASTA format: Uniprot'):
            fasta = FastaEntity(header, sequence, wrong_format, remove_parenthesis)


class TestFastaDataset:
    def _read_all(self, path_list: List[str]) -> str:
        content = ''
        for i in path_list:
            content += Path.read_to_string(i)
        return content

    @pytest.fixture(autouse=True)
    def setup(self):
        # Current Working Directory
        self.CWD = os.path.join(os.getcwd(), 'tests', 'files', 'data', 'fasta')
        # Test Data Folder
        TDF = os.path.join(os.getcwd(), 'tests',  'files', 'data', 'merox')

        fasta_files = Path.list_given_type_files(self.CWD, 'fasta')
        self.fasta_content = self._read_all(fasta_files)
        self.fasta_dataset = FastaDataset(self.fasta_content, 'Custom', False)

        zhrm_folder_path = Path.list_given_type_files(TDF, 'zhrm')
        folder_content = MeroX.load_data(zhrm_folder_path, 'DSBU')
        self.combined_replicas = CrossLink.combine_all(folder_content)

        fas_files = Path.list_given_type_files(self.CWD, 'fas')
        self.fas_content = self._read_all(fas_files)
        self.fas_dataset = FastaDataset(self.fas_content, 'Uniprot')

    def test_positive_constructor(self):
        fasta_dataset = FastaDataset(self.fasta_content, 'Custom')
        assert len(fasta_dataset) == 3

    def test_nonexistent_format_exception_constructor(self):
        with pytest.raises(ValueError, match = 'Wrong FASTA format: Nonexistent'):
            fasta_dataset = FastaDataset(self.fasta_content, 'Nonexistent')

    def test_wrong_format_exception_constructor(self):
        with pytest.raises(ValueError, match = 'Wrong FASTA format: Uniprot'):
            # Custom format is expected
            fasta_dataset = FastaDataset(self.fasta_content, 'Uniprot')

    def test_positive_filter_by_crosslinks(self):
        fasta_dataset = FastaDataset(self.fasta_content, 'Custom', False)
        fasta_dataset.filter_by_crosslinks(self.combined_replicas)
        assert len(fasta_dataset) == 2

    def test_negative_filter_by_crosslinks(self):
        self.combined_replicas.filter_by_score(min_score = sys.maxsize - 1, 
                                          max_score = sys.maxsize)

        self.fasta_dataset.filter_by_crosslinks(self.combined_replicas)
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
        fas_dataset = FastaDataset(self.fas_content, 'Uniprot', True)
        header = r'>sp|P02358|RS6_ECOLI Small ribosomal subunit protein bS6 OS=Escherichia coli strain K12 OX=83333 GN=rpsF PE=1 SV=1'
        fasta_match = fas_dataset.find_gene_by_fasta_header(header)
        expected_match = 'rpsF'
        assert fasta_match == expected_match

    def test_save(self):
        filted_dataset = self.fasta_dataset.filter_by_crosslinks(self.combined_replicas)

        saved_fasta = os.path.join(self.CWD, 'test.fasta')
        if os.path.exists(saved_fasta):
            raise FileExistsError('FastaEntity test file should not exist')

        filted_dataset.save(self.CWD, 'test.fasta')

        if not os.path.exists(saved_fasta):
            raise FileExistsError('FastaEntity test file should exists')

        content = Path.read_to_string(saved_fasta)
        os.remove(saved_fasta)
        saved_fasta_content = FastaDataset(content, 'Custom', False)
        assert str(saved_fasta_content) == str(filted_dataset)


class TestDomainEntity:
    def test_exception_constructor(self):
        input_text = 'wrong format string'
        with pytest.raises(ValueError, match = f'Unknown domain format: {input_text}'):
            domain = DomainEntity(input_text)

    def test_base_color_constructor(self):
        input_text = '>GSTIMP1-cleav,#C0C0C0'
        domain = DomainEntity(input_text)
        assert (
            domain.gene == '>GSTIMP1-cleav' and
            domain.base_color == True and
            domain.color == '#C0C0C0'
        )

    def test_domain_color_constructor(self):
        input_text = '>GSTIMP1-cleav,0,8,#B22222,TAG'
        domain = DomainEntity(input_text)
        assert (
            domain.gene == '>GSTIMP1-cleav' and
            domain.start == 0 and
            domain.end == 8 and
            domain.color == '#B22222' and
            domain.name == 'TAG' and
            domain.base_color == False
        )


class TestDomainDataset:
    def _read_all(self, path_list: List[str]) -> str:
        content = ''
        for i in path_list:
            content += Path.read_to_string(i)
        return content

    @pytest.fixture(autouse=True)
    def setup(self):
        # Current Working Directory
        self.CWD = os.path.join(os.getcwd(), 'tests', 'files', 'data', 'circos_test')
        # Domain Files Directory
        self.DFD = os.path.join(os.getcwd(), 'tests', 'files', 'data', 'dmn')

        domain_path = Path.list_given_type_files(self.DFD, '.dmn')
        domain_content = self._read_all(domain_path)
        self.domains = DomainDataset(domain_content)

    def test_constructor(self):
        assert len(self.domains) == 19

    def test_filter_by_fasta(self):
        fasta_path = os.path.join(os.getcwd(), 'tests', 'files', 'data', 'fasta', 'fasta_1.fasta')
        fasta_content = Path.read_to_string(fasta_path)
        fasta_dataset = FastaDataset(fasta_content, 'Custom')

        len_before = len(self.domains)
        self.domains.filter_by_fasta(fasta_dataset)
        len_after = len(self.domains)
        assert len_before == 19 and len_after == 10
