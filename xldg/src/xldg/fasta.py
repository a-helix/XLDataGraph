from typing import List, Tuple, Iterator
import os
import re

class Fasta_Entity:
    def __init__(self, header: str, sequence: str, fasta_format: str):
        self.raw_header = header.replace('(', '').replace(')', '')  # Merox also removes scopes
        self.raw_sequence = sequence

        if fasta_format == 'Uniprot':
            self.db_id, self.prot_gene = self._split_uniprot_fasta_header(header)
        elif fasta_format == 'Araport11':
            self.db_id, self.prot_gene = self._split_araport11_fasta_header(header)
        elif fasta_format == 'Custom':
            self.db_id, self.prot_gene = self.raw_header, self.raw_header
        else:
            raise ValueError(f'Unknown FASTA format: {fasta_format}')

        #MeroX format of sequence with N-term and C-term as figure brackets
        self.sequence = '{' + self.raw_sequence + '}' 
        self.seq_length = len(self.sequence)

    def _split_uniprot_fasta_header(self, header: str) -> Tuple[str, str]:
        header = header.strip()

        splited_header = header.split('|')
        db_id = splited_header[1]

        prot_gene_match = re.search(r'(GN=[^\s]+)', header)  # Match full 'GN=...'
        prot_gene = prot_gene_match.group(1) if prot_gene_match else ''

        return db_id, prot_gene.replace('GN=', '')

    def _split_araport11_fasta_header(self, header: str) -> Tuple[str, str]:
        splited_header = header.strip().split('|')
        araport11_id = splited_header[0].replace(' ', '').replace('>', '')

        prot_gene = splited_header[1].replace("Symbols: ", "")
        prot_gene = prot_gene.split()
        prot_gene = prot_gene[0].replace(',', '').replace(' ', '')
        return araport11_id, prot_gene

    def __eq__(self, other):
        return (self.raw_header == other.raw_header and  
                self.db_id == other.db_id and 
                self.prot_gene == other.prot_gene and
                self.sequence == other.sequence)

    def __hash__(self):
        return hash((self.raw_header,  
                     self.db_id, 
                     self.prot_gene,
                     self.sequence))

    def __lt__(self, other):
        return self.db_id < other.db_id
    
    def __gt__(self, other):
        return self.db_id > other.db_id


class Fasta_Dataset:
    def __init__(self, fasta_files_paths_list: List[str], fasta_format: str):
        self.fasta_format = fasta_format
        self.entities = self._extract_all_fasta_content_from_folder(fasta_files_paths_list)
        self.size = len(self.entities)
        self._iter_index = 0 
        
    def _extract_all_fasta_content_from_folder(self, fasta_files_path_list: List[str]) -> List['Fasta_Entity']:
        db_entities = []
        
        for file_path in fasta_files_path_list:
            raw_fasta_content = ''
            # Read content of all files
            try:
                with open(file_path, 'r') as file:
                    for line in file:
                        raw_fasta_content += line
            except FileNotFoundError:
                print(f'Fasta_Dataset error: File at {file_path} was not found.')
            except Exception as e:
                print(f'Fasta_Dataset error: {e}')
            
            splited_fasta_content = raw_fasta_content.split('>')
            for fasta in splited_fasta_content:
                splited_fasta = fasta.split('\n')
                
                if len(splited_fasta) > 1:
                    header = '>' + splited_fasta[0]
                    sequence = ''.join(line.strip() for line in splited_fasta[1:])
                    db_entities.append(Fasta_Entity(header, sequence, self.fasta_format))

        sorted_db_entities = sorted(db_entities)
        return sorted_db_entities
    
    def __len__(self):
        return self.size

    def __iter__(self) -> Iterator[Fasta_Entity]:
        self._iter_index = 0
        return self
    
    def __next__(self) -> Fasta_Entity:
        if self._iter_index < self.size:
            entity = self.entities[self._iter_index]
            self._iter_index += 1
            return entity
        else:
            raise StopIteration
        
    def remove_entities_without_merox_xl(self, merox_xls: 'XL_Dataset') -> None:
        filtered_entities = set()
        
        for fasta in self.entities:
            for xl in merox_xls:
                if xl.protein_1 == fasta.raw_header or xl.protein_2 == fasta.raw_header:
                    filtered_entities.add(fasta)
                    break

        # Unifies sector plotting order on a final figure
        self.entities = sorted(list(filtered_entities)) 
        self.size = len(self.entities)

    def find_protein_name_by_header_string(self, header: str) -> str:
        for fasta in self.entities:
            if fasta.raw_header == header:
                return fasta.prot_gene

    def extract_proteins_fasta_enteties_in_XL_Dataset(self, folder_path: str, file_name: str, merox_data: 'XL_Dataset') -> None:
        proteins_in_XL_Dataset = set()

        for xl, _ in merox_data.xls_site_count.items():
            proteins_in_XL_Dataset.add(xl.protein_1)
            proteins_in_XL_Dataset.add(xl.protein_2)

        print(len(proteins_in_XL_Dataset))
        text_output = ''

        print(proteins_in_XL_Dataset)

        for fasta in self.entities:
            print(f'|{fasta.raw_header}|')
            if fasta.raw_header in proteins_in_XL_Dataset:
                text_output += f'{fasta.raw_header}\n'
                text_output += f'{fasta.raw_sequence}\n'
        
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
        
        file_path = os.path.join(folder_path, file_name)
        with open(file_path, "w") as file:
            file.write(text_output)
        
        print(f'Fasta saved to {file_path}')
