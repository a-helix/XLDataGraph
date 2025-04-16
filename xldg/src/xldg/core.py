from dataclasses import dataclass
from typing import List, Tuple, Dict, Set, Iterator, Optional
import os
import sys
import re
import copy
import math
from dataclasses import dataclass
import itertools
import random
import concurrent.futures
import multiprocessing
from functools import partial

import xml.etree.ElementTree as ET
from xml.dom import minidom

import numpy as np
from scipy.spatial import distance
import networkx as nx


class ProteinChainDataset:
    def __init__(self, pcd_data: str):
        self.pcds = {}
        self._assign_pcds(pcd_data)

    def _assign_pcds(self, pcd_data: str) -> None:
        for line in pcd_data.strip().split('\n'):
            splited_line = line.strip().split(',')

            if len(splited_line) < 2:
                raise ValueError(f'Wrong data format in provided input {pcd_data}')

            self.pcds[splited_line[0]] = [item.replace(' ', '') for item in splited_line[1:]]

    def __len__(self):
        return len(self.pcds)

    def __iter__(self):
        return iter(self.pcds.items())

    def __getitem__(self, key):
        return self.pcds[key]

    def __next__(self):
        return next(iter(self.pcds.items()))


class CrossLinkEntity:
    def _remove_text_in_brackets(self, s: str) -> str:
        pattern = r'\(.*?\)'
        cleaned_string = re.sub(pattern, '', s)
        return cleaned_string

    def _initialize_merox_xl(self):
        #Fix for a strange Merox assignment
        self.protein_1 = self._remove_text_in_brackets(self.protein_1).replace('  ', ' ')
        self.from_1 = int(self.from_1)
        self.to_1 = int(self.to_1)
        self.num_site_1 = self.from_1 + int(self.site_1[1:]) - 1

        # Handling N-terminus crosslink
        if self.num_site_1 <= 0:
            self.num_site_1 = 1

        #Fix for a strange Merox assignment
        self.protein_2 = self._remove_text_in_brackets(self.protein_2).replace('  ', ' ')
        self.from_2 = int(self.from_2)
        self.to_2 = int(self.to_2)
        self.num_site_2 =  self.from_2 + int(self.site_2[1:]) - 1

        if self.num_site_2 <= 0:
            self.num_site_2 = 1

        # Handling N-terminus crosslink
        if self.num_site_2 <= 0:
            self.num_site_2 = 1

        self.score = int(self.score)

    def _initialize_predicted_xl(self):
        self.from_1 = int(self.from_1)
        self.to_1 = int(self.to_1)
        self.num_site_1 = int(self.site_1[1:])
        self.from_2 = int(self.from_2)
        self.to_2 = int(self.to_2)
        self.num_site_2 = int(self.site_2[1:])
        self.score = int(self.score)

    def __init__(self, 
                 protein_1: str, 
                 peptide_1: str, 
                 from_1: str, 
                 to_1: str, 
                 site_1: str, 
                 protein_2: str, 
                 peptide_2: str, 
                 from_2: str, 
                 to_2: str, 
                 site_2: str, 
                 score: str, 
                 software: str, 
                 linker: str):

        self.software = software 
        self.linker = linker
        # First peptide data
        self.protein_1 = protein_1 # str
        self.peptide_1 = peptide_1 # str
        self.from_1 = from_1 # int
        self.to_1 = to_1 # int
        self.site_1 = site_1 # str
        self.num_site_1 = None # int
        # Second peptide data
        self.protein_2 = protein_2 # str
        self.peptide_2 = peptide_2 # str
        self.from_2 = from_2 # int
        self.to_2 = to_2 # int
        self.site_2 = site_2 # str
        self.num_site_2 =  None # int
        # Additional info
        self.score = score # int

        if self.software == 'MeroX':
             self._initialize_merox_xl()
        elif self.software == 'Prediction':
            self._initialize_predicted_xl()
        else:
            raise Exception(f'{self.software} is not supported.')

        self.is_interprotein = (self.protein_1 != self.protein_2)
        self.is_homotypical = (self.protein_1 == self.protein_2 and (self.num_site_1 == self.num_site_2 
                                                                     or self.peptide_1 == self.peptide_2))
        if self.software == 'Prediction':
            self.is_homotypical = self.protein_1 == self.protein_2 and self.num_site_1 == self.num_site_2
    
    def __eq__(self, other):
        direct_match = None
        reverse_match = None
        if not isinstance(other, CrossLinkEntity):
            raise ValueError(
                'Cannot compare CrossLinkEntity with non-CrossLinkEntity object')

        if self.software == 'Prediction' or other.software == 'Prediction':
            direct_match = (self.protein_1 == other.protein_1 and
                            self.num_site_1 == other.num_site_1 and
                            self.protein_2 == other.protein_2 and
                            self.num_site_2 == other.num_site_2) 
            reverse_match = (self.protein_1 == other.protein_2 and
                            self.num_site_1 == other.num_site_2 and
                            self.protein_2 == other.protein_1 and
                            self.num_site_2 == other.num_site_1) 
        else:
            direct_match = (self.protein_1 == other.protein_1 and
                            self.peptide_1 == other.peptide_1 and
                            self.site_1 == other.site_1 and
                            self.protein_2 == other.protein_2 and
                            self.peptide_2 == other.peptide_2 and
                            self.site_2 == other.site_2) 
            reverse_match = (self.protein_1 == other.protein_2 and
                            self.peptide_1 == other.peptide_2 and
                            self.site_1 == other.site_2 and
                            self.protein_2 == other.protein_1 and
                            self.peptide_2 == other.peptide_1 and
                            self.site_2 == other.site_1)  
        return direct_match or reverse_match
    
    def __hash__(self):
        if self.protein_1 < self.protein_2 or (self.protein_1 == self.protein_2 and self.num_site_1 <= self.num_site_2):
            return hash((self.protein_1, self.num_site_1, self.protein_2, self.num_site_2))
        else:
            return hash((self.protein_2, self.num_site_2, self.protein_1, self.num_site_1))
    
    def __str__(self):
        return f'{self.protein_1},{self.num_site_1},{self.protein_2},{self.num_site_2},{self.is_interprotein},{self.is_homotypical}'


@dataclass(frozen=True)
class Node:
    '''3D point in space with hashable coordinates'''
    x: float
    y: float
    z: float

    def distance_to(self, other: 'Node') -> float:
        return math.sqrt((self.x - other.x)**2 + 
                         (self.y - other.y)**2 + 
                         (self.z - other.z)**2)

    def to_tuple(self):
        return (self.x, self.y, self.z)


class Atom:
    def __init__(self, number: int, residue: str, type: str, chain: str, x: float, y: float, z: float):
        self.number = number
        self.residue = residue  # Residue identifier (number + insertion code if present)
        self.type = type  # Atom type (e.g., 'N', 'CA')
        self.chain = chain
        self.node = Node(x, y, z)

    def distance_to(self, other: 'Atom') -> float:
        return self.node.distance_to(other.node)

    def __hash__(self):
        return hash((self.number, self.residue, self.type, self.chain))
    
    def __eq__(self, other):
        if not isinstance(other, Atom):
            return False
        return (self.number == other.number and
                self.residue == other.residue and
                self.type == other.type and
                self.chain == other.chain)

    def __lt__(self, other):
        if not isinstance(other, Atom):
            return NotImplemented
        return (self.number, self.residue, self.type, self.chain) < (other.number, other.residue, other.type, other.chain)
    
    def __gt__(self, other):
        if not isinstance(other, Atom):
            return NotImplemented
        return (self.number, self.residue, self.type, self.chain) > (other.number, other.residue, other.type, other.chain)

    def __str__(self):
        return (f'{self.number}\t{self.residue}\t{self.type}\t{self.chain}\t'
                f'{self.node.x:.3f}\t{self.node.y:.3f}\t{self.node.z:.3f}')


class ProteinStructureDataset:
    def __init__(self, file_content: str, format: str):
        self.atoms: List[Atom] = []
        self.primary_nodes = None

        if format == 'pdb':
            self._parse_pdb(file_content)
        elif format == 'cif':
            self._parse_cif(file_content)
        else:
            raise ValueError('Unsupported file format. Only .pdb and .cif formats can be used')

    def _three_to_one_letter_code(self, three_letter_code):
        aa_dict = {
            'ALA': 'A',  # Alanine
            'ARG': 'R',  # Arginine
            'ASN': 'N',  # Asparagine
            'ASP': 'D',  # Aspartic acid
            'CYS': 'C',  # Cysteine
            'GLN': 'Q',  # Glutamine
            'GLU': 'E',  # Glutamic acid
            'GLY': 'G',  # Glycine
            'HIS': 'H',  # Histidine
            'ILE': 'I',  # Isoleucine
            'LEU': 'L',  # Leucine
            'LYS': 'K',  # Lysine
            'MET': 'M',  # Methionine
            'PHE': 'F',  # Phenylalanine
            'PRO': 'P',  # Proline
            'SER': 'S',  # Serine
            'THR': 'T',  # Threonine
            'TRP': 'W',  # Tryptophan
            'TYR': 'Y',  # Tyrosine
            'VAL': 'V',  # Valine
        
            # Non-standard amino acids
            'SEC': 'U',  # Selenocysteine
            'PYL': 'O',  # Pyrrolysine
        
            # Special cases
            'UNK': 'X',  # Unknown
            'XAA': 'X',  # Unspecified/unknown
            'ASX': 'B',  # Asparagine or Aspartic acid
            'GLX': 'Z',  # Glutamine or Glutamic acid
        }
    
        code = three_letter_code.upper()
    
        if code in aa_dict:
            return aa_dict[code]
        else:
            raise ValueError(f'Unknown three-letter amino acid code: {three_letter_code}')

    def _parse_pdb(self, content: str):
        for line in content.split('\n'):
            if line.startswith(('ATOM  ', 'HETATM')):
                try:
                    # Extract fields using fixed column positions
                    atom_number = int(line[6:11].strip())
                    atom_type = line[12:16].strip()
                    residue_name = line[17:20].strip()
                    chain = line[21].strip()
                    residue_number = line[22:26].strip()  # Includes insertion code
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())

                    self.atoms.append(Atom(
                        number = int(residue_number),
                        residue = self._three_to_one_letter_code(residue_name),
                        type = atom_type,
                        chain = chain,
                        x = x,
                        y = y,
                        z = z
                    ))
                except (ValueError, IndexError):
                    continue  

    def _parse_cif(self, content: str):
        atom_data = {}
        current_category = None
        column_names = []
        in_loop = False
    
        # Process file line by line
        lines = content.split('\n')
        i = 0
        while i < len(lines):
            line = lines[i].strip()
        
            # Skip empty lines and comments
            if not line or line.startswith('#'):
                i += 1
                continue
        
            # Check for category
            if line.startswith('_'):
                # Extract category and column name
                parts = line.split('.')
                if len(parts) >= 2:
                    current_category = parts[0].lstrip('_')
                    column_name = parts[1].split()[0]
                
                    # If this is the atom_site category, start collecting column names
                    if current_category == 'atom_site':
                        column_names.append(column_name)
                    
                        # Look ahead to collect all column names
                        j = i + 1
                        while j < len(lines) and lines[j].strip().startswith('_atom_site.'):
                            column_names.append(lines[j].strip().split('.')[1].split()[0])
                            j += 1
                    
                        # Skip the lines we just processed
                        i = j
                        in_loop = True
                        continue
        
            # If we're in the atom_site loop, process the data rows
            elif in_loop and current_category == 'atom_site' and not line.startswith('_') and not line.startswith('#'):
                if line.startswith('loop_'):
                    i += 1
                    continue
                
                # Split the line into fields
                fields = line.split()
            
                # If we've reached the end of the loop or a new category
                if len(fields) == 0 or fields[0].startswith('_') or fields[0].startswith('loop_'):
                    in_loop = False
                    i += 1
                    continue
                
                # Make sure we have enough fields
                if len(fields) >= len(column_names):
                    # Create a dictionary for this row
                    row_data = {column_names[j]: fields[j] for j in range(len(column_names))}
                
                    try:
                        # Extract the required information
                        atom_number = int(row_data.get('id', 0))
                        atom_type = row_data.get('label_atom_id', '').strip()
                        residue_name = row_data.get('label_comp_id', '').strip()
                        chain = row_data.get('auth_asym_id', row_data.get('label_asym_id', '')).strip()
                        residue_number = row_data.get('auth_seq_id', row_data.get('label_seq_id', '')).strip()
                    
                        # Handle insertion code if present
                        pdbx_pdb_ins_code = row_data.get('pdbx_PDB_ins_code', '').strip()
                        if pdbx_pdb_ins_code and pdbx_pdb_ins_code != '?':
                            residue_number += pdbx_pdb_ins_code
                    
                        # Get coordinates
                        x = float(row_data.get('Cartn_x', 0))
                        y = float(row_data.get('Cartn_y', 0))
                        z = float(row_data.get('Cartn_z', 0))
                    
                        # Create and add the atom
                        self.atoms.append(Atom(
                            number = int(residue_number),
                            residue = self._three_to_one_letter_code(residue_name),
                            type = atom_type,
                            chain = chain,
                            x = x,
                            y = y,
                            z = z
                        ))
                    except (ValueError, KeyError) as e:
                        print(f'Error parsing CIF atom entry: {e}')
            i += 1

    def _is_path_distance_in_range(self,
                                   start: int,
                                   goal: int,
                                   min_distance: float,
                                   max_distance: float,
                                   radius: float,
                                   node_multiplier: int
                                   ) -> bool:

        # Step 1: Validate parameters
        start_node = self.primary_nodes[start]
        goal_node = self.primary_nodes[goal]
    
        # Check direct distance first
        direct_distance = start_node.distance_to(goal_node)

        if direct_distance > max_distance:
            return False

        if direct_distance <= max_distance:
            # Check if direct path is valid (no collisions)
            direct_valid = True
            for blocker in self.primary_nodes:
                if self._line_sphere_collision(start_node, goal_node, blocker, radius):
                    direct_valid = False
                    break
        
            if direct_valid and direct_distance >= min_distance:
                return True
        

        # Step 2: Generate helper nodes focused around start and end points
        num_helper_nodes = len(self.primary_nodes) * node_multiplier
        helper_nodes = self._generate_prioritized_helpers(start_node, goal_node, radius, num_helper_nodes, max_distance)
        helper_nodes += [start_node, goal_node]

        # Step 3: Build optimized visibility graph
        G = self._build_optimized_visibility_graph(helper_nodes, radius, max_distance, start_node, goal_node)

        # Step 4: Pathfinding
        try:
            path_length = nx.dijkstra_path_length(G, start_node, goal_node, weight='weight')
            # Check if path length is within range
            return min_distance <= path_length <= max_distance
        except nx.NetworkXNoPath:
            # No valid path exists
            return False

    def _is_point_valid(self, point: Node, obstacles: list[Node], radius: float) -> bool:
        for center in obstacles:
            if point.distance_to(center) < radius:
                return False
        return True

    def _generate_prioritized_helpers(self, 
                                      start_node: Node, 
                                      goal_node: Node, 
                                      radius: float, 
                                      max_samples: int, 
                                      max_distance: float
                                      ) -> list[Node]:
        '''
        Generate helper nodes prioritizing the areas around start and end points
        with compact placement to maximize coverage while minimizing node count.
        '''
        # Calculate midpoint between start and goal
        midpoint = Node(
            (start_node.x + goal_node.x) / 2,
            (start_node.y + goal_node.y) / 2,
            (start_node.z + goal_node.z) / 2
        )
    
        # Calculate direct distance between start and goal
        direct_distance = start_node.distance_to(goal_node)
    
        # Find obstacles within the relevant search area
        relevant_obstacles = []
        extended_search = max_distance * 1.2  # Add a buffer for potential detours
    
        for obs in self.primary_nodes:
            # Check if obstacle is near the path area
            if (obs.distance_to(start_node) <= extended_search or 
                obs.distance_to(goal_node) <= extended_search or
                obs.distance_to(midpoint) <= extended_search):
                relevant_obstacles.append(obs)
    
        # Auto-generate number of samples based on problem constraints
        path_complexity = min(len(relevant_obstacles), 50)
        distance_factor = direct_distance / radius
    
        auto_samples = min(
            max_samples,
            int(10 + path_complexity * 0.5 + distance_factor * 2)
        )
    
        if max_distance > direct_distance * 2:
            auto_samples = min(max_samples, int(auto_samples * 1.5))
    
        auto_samples = max(20, auto_samples)
        num_samples = auto_samples
    
        # Set up bounds for sampling
        coords = np.array([n.to_tuple() for n in relevant_obstacles])
        min_bounds = coords.min(axis=0) - radius * 2
        max_bounds = coords.max(axis=0) + radius * 2
    
        # Define regions for compact sampling
        start_region_radius = min(max_distance * 0.4, direct_distance * 0.6)
        goal_region_radius = min(max_distance * 0.4, direct_distance * 0.6)
        middle_region_radius = min(max_distance * 0.6, direct_distance)
    
        # Allocate samples to regions
        start_region_target = num_samples * 0.35
        goal_region_target = num_samples * 0.35
        middle_region_target = num_samples * 0.3
    
        # Compact placement approach
        # 1. Use grid-based sampling with spacing close to 2*radius
        # 2. Perturb grid points slightly to maximize coverage
        # 3. Accept only valid points
    
        helpers = []
        attempts = 0
        start_region_nodes = 0
        goal_region_nodes = 0
        middle_region_nodes = 0
    
        # Grid spacing for compact placement (slightly larger than 2*radius)
        grid_spacing = radius * 2.2  # Small buffer to ensure minimum separation
    
        # Storage for placed points to facilitate compact placement
        placed_points = []
        max_attempts = num_samples * 15
    
        while (len(helpers) < num_samples and attempts < max_attempts):
            # Determine which region to sample from
            if start_region_nodes < start_region_target:
                center = start_node
                sampling_radius = start_region_radius
                region = 'start'
            elif goal_region_nodes < goal_region_target:
                center = goal_node
                sampling_radius = goal_region_radius
                region = 'goal'
            else:
                center = midpoint
                sampling_radius = middle_region_radius
                region = 'middle'
        
            # Try compact placement - first check if we have existing points
            if placed_points:
                reference_point = placed_points[attempts % len(placed_points)]

                golden_ratio = (1 + 5 ** 0.5) / 2
                theta = (attempts * golden_ratio) % (2 * np.pi)
                phi = (attempts * golden_ratio * 0.5) % np.pi
            
                # Convert to cartesian coordinates on unit sphere
                dx = np.sin(phi) * np.cos(theta)
                dy = np.sin(phi) * np.sin(theta)
                dz = np.cos(phi)
            
                # Scale to grid spacing and offset
                x = reference_point.x + dx * grid_spacing
                y = reference_point.y + dy * grid_spacing
                z = reference_point.z + dz * grid_spacing
            
                # Ensure within region and bounds
                dist_to_center = np.sqrt((x - center.x)**2 + (y - center.y)**2 + (z - center.z)**2)
                if dist_to_center > sampling_radius:
                    # If outside region radius, place on region boundary in that direction
                    scale_factor = sampling_radius / dist_to_center
                    x = center.x + (x - center.x) * scale_factor
                    y = center.y + (y - center.y) * scale_factor
                    z = center.z + (z - center.z) * scale_factor
            
                # Ensure within global bounds
                x = max(min_bounds[0], min(max_bounds[0], x))
                y = max(min_bounds[1], min(max_bounds[1], y))
                z = max(min_bounds[2], min(max_bounds[2], z))
            
                point = Node(x, y, z)
            else:
                # No existing points yet, start with biased sampling near center
                point = self._sample_around_point(center, sampling_radius * 0.5, min_bounds, max_bounds, attempts)
        
            # Validate point
            valid = self._is_point_valid(point, relevant_obstacles, radius)
        
            # Also check minimum distance from other helper nodes for compact placement
            if valid:
                for existing in helpers:
                    if point.distance_to(existing) < radius * 2:
                        valid = False
                        break
        
            if valid:
                helpers.append(point)
                placed_points.append(point)
            
                # Update region counters
                if region == 'start':
                    start_region_nodes += 1
                elif region == 'goal':
                    goal_region_nodes += 1
                else:
                    middle_region_nodes += 1
        
            attempts += 1
        
            # Adaptive adjustment
            if attempts > max_attempts * 0.7 and len(helpers) < num_samples * 0.5:
                # If placement is difficult, gradually reduce spacing requirements
                grid_spacing = max(radius * 1.8, grid_spacing * 0.95)
            
                # Also reduce targets
                num_samples = max(len(helpers) + 5, int(num_samples * 0.8))
                start_region_target = num_samples * 0.35
                goal_region_target = num_samples * 0.35
                middle_region_target = num_samples * 0.3
    
        return helpers

    def _sample_around_point(self, center, sampling_radius, min_bounds, max_bounds, attempt):
        golden_ratio = (1 + 5 ** 0.5) / 2
        theta = (attempt * golden_ratio) % (2 * np.pi)
        phi = (attempt * golden_ratio * 0.5) % np.pi

        dx = np.sin(phi) * np.cos(theta)
        dy = np.sin(phi) * np.sin(theta)
        dz = np.cos(phi)

        x = center.x + dx * sampling_radius
        y = center.y + dy * sampling_radius
        z = center.z + dz * sampling_radius

        x = max(min_bounds[0], min(max_bounds[0], x))
        y = max(min_bounds[1], min(max_bounds[1], y))
        z = max(min_bounds[2], min(max_bounds[2], z))

        return Node(x, y, z)

    def _line_sphere_collision(self, p1: Node, p2: Node, sphere_center: Node, sphere_radius: float) -> bool:
        p1 = np.array(p1.to_tuple())
        p2 = np.array(p2.to_tuple())
        c = np.array(sphere_center.to_tuple())
    
        # Line segment direction vector
        d = p2 - p1
    
        # Vector from sphere center to first point
        f = p1 - c
    
        # Quadratic equation coefficients: at² + bt + c = 0
        a = np.dot(d, d)  # Squared length of line segment
        b = 2 * np.dot(f, d)
        c_val = np.dot(f, f) - sphere_radius ** 2
    
        # Calculate discriminant
        discriminant = b ** 2 - 4 * a * c_val
    
        # No intersection if discriminant is negative
        if discriminant < 0:
            return False
        
        # Calculate intersection points
        discriminant = np.sqrt(discriminant)
        t1 = (-b - discriminant) / (2 * a)
        t2 = (-b + discriminant) / (2 * a)
    
        # Line segment intersects sphere if at least one intersection point
        # lies on the segment (0 <= t <= 1)
        return (0 <= t1 <= 1) or (0 <= t2 <= 1)

    def _build_optimized_visibility_graph(self, 
                                          helper_nodes: list[Node], 
                                          radius: float,
                                          max_distance: float,
                                          start_node: Node,
                                          goal_node: Node
                                          ):
        '''Build visibility graph optimized to focus on nodes that might contribute to a valid path'''
        G = nx.Graph()
    
        # Add all nodes to the graph
        for node in helper_nodes:
            G.add_node(node)
    
        # Calculate relevance of primary nodes (obstacles)
        relevant_obstacles = []
        extended_search = max_distance * 1.5  # Add buffer for potential detours
    
        for obs in self.primary_nodes:
            # Check if obstacle is near the path area
            if (obs.distance_to(start_node) <= extended_search and 
                obs.distance_to(goal_node) <= extended_search):
                relevant_obstacles.append(obs)
    
        # Sort helper nodes by proximity to start and goal
        sorted_helpers = sorted(
            helper_nodes,
            key=lambda node: min(node.distance_to(start_node), node.distance_to(goal_node))
        )
    
        # Connect start and goal to nearby helpers first
        connect_radius = max_distance * 0.5  # Connect to nodes within this radius
    
        # Connect start node to nearby helpers
        for helper in sorted_helpers:
            if helper != start_node and helper != goal_node:
                dist = start_node.distance_to(helper)
                if dist <= connect_radius:
                    collision = False
                    for blocker in relevant_obstacles:
                        if blocker != start_node and self._line_sphere_collision(start_node, helper, blocker, radius):
                            collision = True
                            break
                    if not collision:
                        G.add_edge(start_node, helper, weight=dist)
    
        # Connect goal node to nearby helpers
        for helper in sorted_helpers:
            if helper != start_node and helper != goal_node:
                dist = goal_node.distance_to(helper)
                if dist <= connect_radius:
                    collision = False
                    for blocker in relevant_obstacles:
                        if blocker != goal_node and self._line_sphere_collision(goal_node, helper, blocker, radius):
                            collision = True
                            break
                    if not collision:
                        G.add_edge(goal_node, helper, weight=dist)
    
        # Connect other helper nodes, focusing on those closer to start/goal first
        for i, p1 in enumerate(sorted_helpers):
            # Skip if already connected to both start and goal
            if p1 != start_node and p1 != goal_node:
                # Only check connections with nodes that are within a reasonable distance
                for j in range(i + 1, len(sorted_helpers)):
                    p2 = sorted_helpers[j]
                    # Quick distance check before expensive collision check
                    dist = p1.distance_to(p2)
                
                    # Only connect if they're close enough (but not too close)
                    if dist <= connect_radius * 1.5:
                        collision = False
                        for blocker in relevant_obstacles:
                            if self._line_sphere_collision(p1, p2, blocker, radius):
                                collision = True
                                break
                        if not collision:
                            G.add_edge(p1, p2, weight=dist)
        return G

    def _check_path_distance_mp(self, args):
        # Unpack the arguments
        psd_instance, idx1, idx2, min_length, max_length, radius, node_multiplier = args
    
        # Call the instance method
        result = psd_instance._is_path_distance_in_range(
            idx1, idx2, min_length, max_length, radius, node_multiplier
        )
        return (idx1, idx2, result)

    def predict_crosslinks(self, 
                           pcd: ProteinChainDataset, 
                           residues_1: str, 
                           residues_2: str, 
                           min_length: float = 1.0, 
                           max_length: float = sys.maxsize, 
                           linker: Optional[str] = None,
                           atom_type: str = 'CA',
                           direct_path = True,
                           radius: float = 1.925, # Half of CA-CA distance in peptide bond which is typicall between 3.75 to 3.85
                           node_multiplier: int = 100,
                           num_processes: int = 1 
                           ) -> 'CrossLinkDataset':

        filtered_atoms = list(filter(lambda atom: atom.type == atom_type, self.atoms))

        if len(filtered_atoms) == 0:
            raise ValueError(f'Invalid atom type {atom_type}')

        if self.primary_nodes is None:
            self.primary_nodes = [atom.node for atom in self.atoms if atom.type == atom_type]

        if node_multiplier < 1:
            raise ValueError(
                f'Node multiplier must be greater than 0, got {node_multiplier}')

        if min_length <= 0:
            raise ValueError(
                f'Min length {min_length} must be greater than 0')

        if max_length < min_length:
            raise ValueError(
                f'Max length {max_length} must be greater than min length {min_length}')
    
        # Set reasonable process count based on system
        if num_processes <= 0:
            num_processes = max(1, multiprocessing.cpu_count())

        def _get_all_crosslink_candidates(residues: str) -> List['Atom']:
            # Same as original implementation
            residue_chars = list(residues)
            atoms = []
            for _, chains in pcd:
                for chain in chains:
                    for res in residue_chars:
                        if res == '{':
                            # Find N-terminus atom
                            for atom in filtered_atoms:
                                if atom.chain == chain and atom.type == atom_type and atom.number == 1:
                                    atoms.append(atom)
                                    break
                        elif res == '}':
                            # Find C-terminus atom
                            last_atom = None
                            counter = 0
                            for atom in filtered_atoms:
                                if atom.chain == chain and atom.type == atom_type and atom.number > counter:
                                    counter = atom.number
                                    last_atom = atom
                            if last_atom:  # Only append if we found an atom
                                atoms.append(last_atom)
                        else:
                            # Find all atoms of specific residue type
                            for atom in filtered_atoms:
                                if atom.chain == chain and atom.type == atom_type and atom.residue == res:
                                    atoms.append(atom)
            return atoms

        atoms_1 = _get_all_crosslink_candidates(residues_1)
        atoms_2 = _get_all_crosslink_candidates(residues_2)

        crosslink_pairs = []  # Using list instead of set to preserve order

        if direct_path:
            # Direct path calculation (not parallelized as it's already fast)
            for a1 in atoms_1:
                for a2 in atoms_2:
                    # Avoid self-crosslinks (same atom)
                    if a1 is a2:
                        continue
                    if min_length <= a1.distance_to(a2) <= max_length:
                        crosslink_pairs.append((a1, a2))
        else:
            # Path distance calculation (multi-process version)
            # Generate all pairs to evaluate and convert to indices
            pairs_to_check = []
            for a1 in atoms_1:
                for a2 in atoms_2:
                    # Avoid self-crosslinks (same atom)
                    if a1 is not a2:
                        idx1 = filtered_atoms.index(a1)
                        idx2 = filtered_atoms.index(a2)
                        # Create input tuple with all parameters needed
                        pairs_to_check.append((self, idx1, idx2, min_length, max_length, radius, node_multiplier))
        
            # Process in parallel using Pool
            pool = multiprocessing.Pool(processes=num_processes)
            try:
                # Use map to parallelize
                results = pool.map(self._check_path_distance_mp, pairs_to_check, chunksize=10)
            
                # Filter positive results and convert back to atom pairs
                for idx1, idx2, is_in_range in results:
                    if is_in_range:
                        crosslink_pairs.append((filtered_atoms[idx1], filtered_atoms[idx2]))
            finally:
                pool.close()
                pool.join()

        def _get_protein_by_chain(chain: str) -> str:
            for protein, chains in pcd:
                for ch in chains:
                    if ch == chain:
                        return protein
            raise ValueError(f'Undefined protein chain {chain}')

        crosslinks = []
        for a1, a2 in crosslink_pairs:
            crosslinks.append(CrossLinkEntity(
                _get_protein_by_chain(a1.chain),
                '',
                '-1',
                '-1',
                a1.residue + str(a1.number),
                _get_protein_by_chain(a2.chain),
                '',
                '-1',
                '-1',
                a2.residue + str(a2.number),
                '-1',
                'Prediction',
                linker
            ))
        result = CrossLinkDataset(crosslinks)
        return result.blank_replica_counter()
        
    def __iter__(self):
        return iter(self.atoms)

    def __len__(self):
        return len(self.atoms)

class CrossLinkDataset:
    def __init__(self, xls: List['CrossLinkEntity']):
        self.xls = xls
        self._remove_invalid_xls()

        self.xls_site_count = self._quantify_elements(self.xls)

    def __iter__(self):
        self._index = 0
        return self

    def __next__(self):
        if self._index < len(self.xls):
            result = self.xls[self._index]
            self._index += 1
            return result
        else:
            raise StopIteration

    def __add__(self, other):
        if not isinstance(other, CrossLinkDataset):
            return NotImplemented  # Return if other is not a CrossLinkDataset

        combined_xls = self.xls + other.xls
        combined_xls_site_count = self.xls_site_count

        for site, count in other.xls_site_count.items():
            if site not in combined_xls_site_count:
                combined_xls_site_count[site] = count 
            else:
                combined_xls_site_count[site] += count

        final = CrossLinkDataset(combined_xls)
        final.xls_site_count = combined_xls_site_count

        return final

    def __getitem__(self, index):
        return self.xls[index]

    def __len__(self):
        return len(self.xls)

    def filter_by_score(self, min_score: int = 0, max_score: int = sys.maxsize) -> 'CrossLinkDataset':
        filtered_list = []

        for xl in self.xls:
            if xl.score >= min_score and xl.score <= max_score:
                filtered_list.append(xl)
        
        unique_filtered_list = set(filtered_list)
        filterered_xls_site_count = {}

        for xl, count in self.xls_site_count.items():
            if xl in unique_filtered_list:
                filterered_xls_site_count[xl] = count
                    
        self.xls = filtered_list
        self.xls_site_count = filterered_xls_site_count
        return self

    def filter_by_replica(self, min_rep: int = 1, max_rep: int = sys.maxsize) -> 'CrossLinkDataset':
        if min_rep > max_rep:
            raise ValueError('CircosConfig max_rep cannot be smaller than min_rep')

        filtered_xls_site_count = {}
        filtered_xls = []

        for xl1, count in self.xls_site_count.items():
            if count >= min_rep and count <= max_rep:
                filtered_xls_site_count[xl1] = count
                for xl2 in self.xls:
                    if xl2 == xl1:
                        filtered_xls.append(xl2)

        self.xls_site_count = filtered_xls_site_count
        self.xls = filtered_xls
        return self

    def remove_interprotein_crosslinks(self) -> 'CrossLinkDataset':
        filtered_xls = []
        for xl in self.xls:
            if xl.is_homotypical:
                filtered_xls.append(xl)
                continue
            if xl.is_interprotein is False:
                filtered_xls.append(xl)

        self._update_xls_data(filtered_xls)
        return self

    def remove_intraprotein_crosslinks(self) -> 'CrossLinkDataset':
        filtered_xls = []
        for xl in self.xls:
            if xl.is_homotypical:
                filtered_xls.append(xl)
                continue
            if xl.is_interprotein is True:
                filtered_xls.append(xl)

        self._update_xls_data(filtered_xls)
        return self

    def remove_homotypic_crosslinks(self) -> 'CrossLinkDataset':
        filtered_xls = []
        for xl in self.xls:
            if xl.is_homotypical is False:
                filtered_xls.append(xl)

        self._update_xls_data(filtered_xls)
        return self

    def _update_xls_data(self, xls: List['CrossLinkEntity']) -> None:
        filtered_xls_site_count = {}
        for xl1 in xls:
            for xl2, count in self.xls_site_count.items():
                if xl1 == xl2:
                    filtered_xls_site_count[xl2] = count

        self.xls = xls
        self.xls_site_count = filtered_xls_site_count

    def blank_replica_counter(self)  -> 'CrossLinkDataset':
        for key in self.xls_site_count.keys():
            self.xls_site_count[key] = 1  
        return self
        
    def _quantify_elements(self, elements: List['CrossLinkEntity']) -> Dict['CrossLinkEntity', int]:
        element_counts = {}
        for element in elements:
            if element not in element_counts:
                element_counts[element] = 1 
            else:
                element_counts[element] += 1

        return element_counts

    def _remove_invalid_xls(self) -> None:
        buffer = []
        for xl in self.xls:
            if xl.software == 'MeroX':
                # Ignore MeroX decoy matches
                if 'DEC_' in xl.protein_1 or 'DEC_' in xl.protein_2:
                    continue
                # Ignore MeroX dead-end matches
                if 'x0' in xl.site_1 or 'x0' in xl.site_2:
                    continue
                
                if 'H2O' in xl.protein_1 or 'H2O' in xl.protein_2:
                    continue

            buffer.append(xl)
        self.xls = buffer

    def save_crosslink_counter(self, folder_path: str, file_name: str, separator: str = '\t'):
        file = os.path.join(folder_path, file_name)
        os.makedirs(folder_path, exist_ok=True)

        str_file = file
        header = (
            f'software{separator}protein_1{separator}peptide_1{separator}from_1{separator}'
            f'to_1{separator}site_1{separator}protein_2{separator}peptide_2{separator}'
            f'from_2{separator}to_2{separator}site_2{separator}interprotein{separator}'
            f'homotypical{separator}replicas\n'
        )

        with open(file, 'w') as file:
            file.write(header)
            for xl, frequency in self.xls_site_count.items():
                file.write(
                    f'{xl.software}{separator}{xl.protein_1}{separator}{xl.peptide_1}{separator}'
                    f'{xl.from_1}{separator}{xl.to_1}{separator}{xl.site_1}{separator}'
                    f'{xl.protein_2}{separator}{xl.peptide_2}{separator}{xl.from_2}{separator}'
                    f'{xl.to_2}{separator}{xl.site_2}{separator}{xl.is_interprotein}{separator}'
                    f'{xl.is_homotypical}{separator}{frequency}\n'
                )

    def save_crosslinks(self, folder_path: str, file_name: str, separator: str = '\t'):
        file = os.path.join(folder_path, file_name)
        os.makedirs(folder_path, exist_ok=True)

        str_file = file
        header = (
            f'software{separator}protein_1{separator}peptide_1{separator}from_1{separator}'
            f'to_1{separator}site_1{separator}protein_2{separator}peptide_2{separator}'
            f'from_2{separator}to_2{separator}site_2{separator}interprotein{separator}'
            f'homotypical{separator}score\n'
        )

        with open(file, 'w') as file:
            file.write(header)
            for xl in self.xls:
                file.write(
                    f'{xl.software}{separator}{xl.protein_1}{separator}{xl.peptide_1}{separator}'
                    f'{xl.from_1}{separator}{xl.to_1}{separator}{xl.site_1}{separator}'
                    f'{xl.protein_2}{separator}{xl.peptide_2}{separator}{xl.from_2}{separator}'
                    f'{xl.to_2}{separator}{xl.site_2}{separator}{xl.is_interprotein}{separator}'
                    f'{xl.is_homotypical}{separator}{xl.score}\n'
                )

    def export_for_chimerax(self, 
                            pcd: ProteinChainDataset, 
                            folder_path: str, 
                            file_name: str,  
                            diameter: float = 0.2, 
                            dashes: int = 1,
                            color_valid_distance: str = '#48cae4', # Sky Blue
                            color_invalid_outsider: str = '#d62828', # Red
                            protein_structure: ProteinStructureDataset = None,
                            min_distance: float = 0,
                            max_distance: float = sys.maxsize,
                            atom_type: str = 'CA'
                            ) -> None:

        os.makedirs(folder_path, exist_ok=True)
    
        # Only filter protein_structure if it's not None
        filtered_protein_structure = []
        if protein_structure is not None:
            filtered_protein_structure = list(filter(lambda atom: atom.type == atom_type, protein_structure.atoms))
            if len(filtered_protein_structure) == 0:
                raise ValueError(f'Invalid atom type {atom_type}')

        def _write_to_pb_file(buffer: List[str], filename: str):
            buffer = '\n'.join(buffer)
            if buffer:
                file_path = os.path.join(folder_path, filename)
                with open(file_path, 'w') as file:
                    file.write(f'; dashes = {dashes}\n; radius = {diameter}\n{buffer}')

        def _clasify_crosslink(crosslink: CrossLinkEntity, valid_buffer: List[str], outliar_buffer: List[str], chain1: str, chain2: str) -> None:
            # Create lookup dictionaries for faster access if protein_structure exists
            if filtered_protein_structure:
                # Convert site numbers to strings for comparison if needed
                site1_str = str(crosslink.num_site_1)
                site2_str = str(crosslink.num_site_2)
            
                # Find matching atoms for both sites
                matching_atoms1 = [atom for atom in filtered_protein_structure if str(atom.number) == site1_str and atom.chain == chain1]
                matching_atoms2 = [atom for atom in filtered_protein_structure if str(atom.number) == site2_str and atom.chain == chain2]

                if len(matching_atoms1) != 1 or len(matching_atoms1) != 1:
                    raise ValueError(f'Wrong chain assignment in ProteinChain file')

                distance = matching_atoms1[0].distance_to(matching_atoms2[0])
                if distance >= min_distance and distance <= max_distance:
                    valid_buffer.append(f'/{chain1}:{crosslink.num_site_1}@CA\t/{chain2}:{crosslink.num_site_2}@CA\t{color_valid_distance}')
                else:
                    outliar_buffer.append(f'/{chain1}:{crosslink.num_site_1}@CA\t/{chain2}:{crosslink.num_site_2}@CA\t{color_invalid_outsider}')

            else:
                # No protein structure to validate against, consider all valid
                valid_buffer.append(f'/{chain1}:{crosslink.num_site_1}@CA\t/{chain2}:{crosslink.num_site_2}@CA\t{color_valid_distance}')

        xl_frequencies: Set[int] = set(self.xls_site_count.values())
        for xl_frequency in xl_frequencies:
            buffer_homotypical_xl = []
            buffer_heterotypical_intra_xl = []
            buffer_heterotypical_inter_xl = []

            outliers_buffer_homotypical_xl = []
            outliers_buffer_heterotypical_intra_xl = []
            outliers_buffer_heterotypical_inter_xl = []
            distance = 0

            for key, value in self.xls_site_count.items():
                if value != xl_frequency:
                    continue

                crosslink = self._validate_terminus_sites(key)
                chains = pcd[crosslink.protein_1] 
                if crosslink.is_homotypical:
                    for c1 in chains:
                        for c2 in chains:
                            if c1 != c2:  # Ensure unique chain pairing
                                _clasify_crosslink(crosslink, buffer_homotypical_xl, outliers_buffer_homotypical_xl, c1, c2)

                elif crosslink.is_interprotein:
                    chain1 = pcd[crosslink.protein_1]
                    chain2 = pcd[crosslink.protein_2]

                    for c1 in chain1:
                        for c2 in chain2:
                            _clasify_crosslink(crosslink, buffer_heterotypical_inter_xl, outliers_buffer_heterotypical_inter_xl, c1, c2)

                else:  # Intraprotein heterotypical
                    for c1 in chains:
                        for c2 in chains:
                            _clasify_crosslink(crosslink, buffer_heterotypical_intra_xl, outliers_buffer_heterotypical_intra_xl, c1, c2)

        # Writing to files
        if protein_structure:
            _write_to_pb_file(outliers_buffer_heterotypical_inter_xl, f'outliers_{file_name}_interprotein_xl_{xl_frequency}_rep.pb')
            _write_to_pb_file(outliers_buffer_heterotypical_intra_xl, f'outliers_{file_name}_intraprotein_xl_{xl_frequency}_rep.pb')
            _write_to_pb_file(outliers_buffer_homotypical_xl, f'outliers_{file_name}_homotypical_xl_{xl_frequency}_rep.pb')

        _write_to_pb_file(buffer_heterotypical_inter_xl, f'{file_name}_interprotein_xl_{xl_frequency}_rep.pb')
        _write_to_pb_file(buffer_heterotypical_intra_xl, f'{file_name}_intraprotein_xl_{xl_frequency}_rep.pb')
        _write_to_pb_file(buffer_homotypical_xl, f'{file_name}_homotypical_xl_{xl_frequency}_rep.pb')

    def export_ppis_for_gephi(
        self, 
        pcd: ProteinChainDataset, 
        folder_path: str, 
        file_name: str
    ) -> None:

        save_path = self._validate_gephi_format(folder_path, file_name)
        node_buffer = dict()
        edge_buffer = dict()
    
        for key, value in pcd:
            for i in value:
                node_buffer[i] = key
    
        for xl in self.xls_site_count:
            for chain1, id1 in node_buffer.items():
                for chain2, id2 in node_buffer.items():
                    if chain1 == chain2:
                        continue
                    xl = self._validate_terminus_sites(xl)
                    if xl.protein_1 == id1 and xl.protein_2 == id2:
                        pair = tuple(sorted((chain1, chain2)))
                        edge_buffer[pair] = edge_buffer.get(pair, 0) + 1
    
        self._create_gexf(save_path, node_buffer, edge_buffer)

    def export_aais_for_gephi(self, 
                              pcd: ProteinChainDataset, 
                              folder_path: str, 
                              file_name: str 
                              ) -> None:

        save_path = self._validate_gephi_format(folder_path, file_name)
        node_buffer = dict()
        edge_buffer = dict()
        size_buffer = dict()

        # Calculate min/max sites for each protein
        for key, _ in pcd:
            min_site = sys.maxsize
            max_site = 0

            for xl in self.xls_site_count:
                validated_xl = self._validate_terminus_sites(xl)
                if validated_xl.protein_1 == key:
                    min_site = min(min_site, validated_xl.num_site_1)
                    max_site = max(max_site, validated_xl.num_site_1)
                if validated_xl.protein_2 == key:
                    min_site = min(min_site, validated_xl.num_site_2)
                    max_site = max(max_site, validated_xl.num_site_2)

            size_buffer[key] = (min_site, max_site)

        for key, value in pcd:
            for chain in value:
                first = size_buffer[key][0]
                last = size_buffer[key][1] + 1
                for i in range(first, last):
                    node_id = f'{chain}{i}_{key}'
                    node_label = f'{chain}_{i}'
                    node_buffer[node_id] = node_label

                for i in range(first, last - 1):
                    node_id = f'{chain}{i}_{key}'
                    next_node = f'{chain}{i+1}_{key}'
                    edge_buffer[(node_id, next_node)] = 1

        # Map chains to their protein IDs
        chain_buffer = {chain: key for key, chains in pcd for chain in chains}

        # Process CrossLinkEntitys between chains
        for xl in self.xls_site_count:
            validated_xl = self._validate_terminus_sites(xl)
            for chain1, id1 in chain_buffer.items():
                for chain2, id2 in chain_buffer.items():

                    if (validated_xl.protein_1 == id1 and validated_xl.protein_2 == id2):
                        first_id = f'{chain1}{validated_xl.num_site_1}_{id1}'
                        second_id = f'{chain2}{validated_xl.num_site_2}_{id2}'
                        pair = tuple(sorted((first_id, second_id)))
                        edge_buffer[pair] = edge_buffer.get(pair, 0) + 1

        self._create_gexf(save_path, node_buffer, edge_buffer)

    def _validate_terminus_sites(self, CrossLinkEntity: 'CrossLinkEntity') -> 'CrossLinkEntity':
        buffer = copy.deepcopy(CrossLinkEntity)
        if '{' in buffer.site_1:
            buffer.num_site_1 += 1 
        if '}' in buffer.site_1:
            buffer.num_site_1 -= 1 
        if '{' in buffer.site_2:
            buffer.num_site_2 += 1 
        if '}' in buffer.site_2:
            buffer.num_site_2 -= 1 
        return buffer

    def _validate_gephi_format(self, folder_name: str, file_name: str) -> str:
        file_extension = os.path.splitext(file_name)[1]
        lower_extension = file_extension.lower()
        if lower_extension != '.gexf':
            raise ValueError(f'Wrong data format is provided in {file_name}. Only ".gexf" format is supported')
        os.makedirs(folder_name, exist_ok=True)
        return os.path.join(folder_name, file_name)

    def _create_gexf(self, file_path: str, node_buffer, edge_buffer) -> None:
        gexf = ET.Element('gexf', version='1.3', xmlns='http://www.gexf.net/1.3')
        graph = ET.SubElement(gexf, 'graph', mode='static', defaultedgetype='undirected')
    
        nodes = ET.SubElement(graph, 'nodes')
        for node_id, label in node_buffer.items():
            ET.SubElement(nodes, 'node', {'id': str(node_id), 'label': str(label)})
    
        edges = ET.SubElement(graph, 'edges')
        for edge_id, ((source, target), weight) in enumerate(edge_buffer.items()):
            ET.SubElement(edges, 'edge', {
                'id': str(edge_id),
                'source': str(source),
                'target': str(target),
                'weight': str(weight)
            })
    
        rough_string = ET.tostring(gexf, encoding='utf-8')
        reparsed = minidom.parseString(rough_string)
        xml_content = reparsed.toprettyxml(indent='  ')
        xml_content = xml_content.replace('<?xml version="1.0" ?>',
                            '<?xml version="1.0" encoding="UTF-8"?>')

        with open(file_path, 'w', encoding='utf-8') as f:
            f.write(xml_content)

    @classmethod
    def unique_elements(cls, dataset1: 'CrossLinkDataset', dataset2: 'CrossLinkDataset') -> Tuple['CrossLinkDataset', 'CrossLinkDataset']:
        count1 = dataset1.xls_site_count
        count2 = dataset1.xls_site_count
        
        set1 = set(dataset1.xls)
        set2 = set(dataset2.xls)
        
        unique_to_dataset1 = [xl for xl in dataset1.xls if xl not in set2]
        unique_to_dataset2 = [xl for xl in dataset2.xls if xl not in set1]
        
        # Create datasets from unique elements
        unique_dataset1 = cls(unique_to_dataset1)
        unique_dataset2 = cls(unique_to_dataset2)

        def _filter_xls_site_count(xls_site_count: Dict['CrossLinkEntity', int], list_of_xls: List['CrossLinkEntity']) -> 'CrossLinkDataset':
            return {k: xls_site_count[k] for k in xls_site_count if k in list_of_xls}

        # Set the xls_site_count for the unique datasets
        unique_dataset1.xls_site_count = _filter_xls_site_count(count1, unique_to_dataset1)
        unique_dataset2.xls_site_count = _filter_xls_site_count(count2, unique_to_dataset2)
        
        return unique_dataset1, unique_dataset2

    @classmethod
    def common_elements(cls, dataset1: 'CrossLinkDataset', dataset2: 'CrossLinkDataset') -> Tuple['CrossLinkDataset', 'CrossLinkDataset']:
        count1 = dataset1.xls_site_count
        count2 = dataset2.xls_site_count
        
        common_elements = set(dataset1.xls) & set(dataset2.xls)
        
        common_list1 = [xl for xl in dataset1.xls if xl in common_elements]
        common_list2 = [xl for xl in dataset2.xls if xl in common_elements]
        
        common_dataset1 = cls(common_list1)
        common_dataset2 = cls(common_list2)

        common_dataset1.xls_site_count = {k: count1[k] for k in count1 if k in common_elements}
        common_dataset2.xls_site_count = {k: count2[k] for k in count2 if k in common_elements}
        
        return common_dataset1, common_dataset2

class FastaEntity:
    def __init__(self, header: str, sequence: str, fasta_format: str, remove_parenthesis: bool = False):
        self.raw_header = header
        self.remove_parenthesis = remove_parenthesis

        if self.remove_parenthesis:
            self.raw_header = header.replace('(', '').replace(')', '')  # Merox also removes scopes

        self.raw_sequence = sequence
        try:
            if fasta_format == 'Uniprot':
                self.db_id, self.prot_gene = self._split_uniprot_fasta_header(header)
            elif fasta_format == 'Araport11':
                self.db_id, self.prot_gene = self._split_araport11_fasta_header(header)
            elif fasta_format == 'Custom':
                self.db_id, self.prot_gene = self.raw_header, self.raw_header
            else:
                raise Exception
        except Exception:
            raise ValueError(f'Wrong FASTA format: {fasta_format}')

                

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

        prot_gene = splited_header[1].replace('Symbols: ', '')
        prot_gene = prot_gene.split()
        prot_gene = prot_gene[0].replace(',', '').replace(' ', '')
        return araport11_id, prot_gene

    def __eq__(self, other):
        return (self.raw_header == other.raw_header and  
                self.db_id == other.db_id and 
                self.prot_gene == other.prot_gene)

    def __hash__(self):
        return hash((self.raw_header,  
                     self.db_id, 
                     self.prot_gene))

    def __lt__(self, other):
        return self.db_id < other.db_id
    
    def __gt__(self, other):
        return self.db_id > other.db_id

    def __str__(self) -> str:
        return f'{self.raw_header}\n{self.raw_sequence}'


class FastaDataset:
    def __init__(self, raw_fasta_content: str, fasta_format: str, remove_parenthesis: bool = False):
        self.fasta_format = fasta_format
        self.remove_parenthesis = remove_parenthesis
        self.entities = self._parse_fasta_content(raw_fasta_content)
        self._iter_index = 0 

    def _parse_fasta_content(self, raw_fasta_content: str) -> List['FastaEntity']:
        db_entities: Set['FastaEntity'] = set()

        if raw_fasta_content:
            db_entities.update(self._process_fasta_entries(raw_fasta_content))

        return sorted(db_entities)

    def _process_fasta_entries(self, raw_fasta_content: str) -> List['FastaEntity']:
        entries = []
        current_header = None
        current_sequence = []

        for line in raw_fasta_content.splitlines():
            line = line.strip()
            if not line:
                continue

            if line.startswith('>'):
                if current_header:
                    header = '>' + current_header
                    if self.remove_parenthesis:
                        header = header.replace(')', '').replace('(', '')

                    entries.append(FastaEntity(header, ''.join(current_sequence), self.fasta_format))
                current_header = line[1:].strip()
                current_sequence = []
            else:
                current_sequence.append(line)

        if current_header:
            header = '>' + current_header
            if self.remove_parenthesis:
                        header = header.replace(')', '').replace('(', '')

            entries.append(FastaEntity(header, ''.join(current_sequence), self.fasta_format))

        return entries
    
    def __len__(self):
        return len(self.entities)

    def __iter__(self) -> Iterator['FastaEntity']:
        self._iter_index = 0
        return self
    
    def __next__(self) -> 'FastaEntity':
        if self._iter_index < len(self.entities):
            entity = self.entities[self._iter_index]
            self._iter_index += 1
            return entity
        else:
            raise StopIteration

    def __str__(self) -> str:
        return '\n'.join(str(fasta) for fasta in self.entities)
        
    def filter_by_crosslinks(self, merox_xls: 'CrossLinkDataset') -> 'FastaDataset':
        filtered_entities = set()
        
        for fasta in self.entities:
            for xl in merox_xls:
                if xl.protein_1 == fasta.raw_header or xl.protein_2 == fasta.raw_header:
                    filtered_entities.add(fasta)
                    break

        # Unifies sector plotting order on a final figure
        self.entities = sorted(list(filtered_entities)) 
        return self

    def find_gene_by_fasta_header(self, header: str) -> str:
        for fasta in self.entities:
            if fasta.raw_header == header:
                return fasta.prot_gene

    def save(self, folder_path: str, file_name: str) -> None:
        text_output = self.__str__()
        
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
        
        file_path = os.path.join(folder_path, file_name)
        with open(file_path, 'w') as file:
            file.write(text_output)


class DomainEntity:
    def __init__(self, input: str):
        self.splited_data = input.split(',')
        if len(self.splited_data) == 5:
            self.gene = self.splited_data[0].replace(' ', '')
            self.start = int(self.splited_data[1].replace(' ', ''))
            self.end = int(self.splited_data[2].replace(' ', ''))
            self.color = self.splited_data[3].replace(' ', '')
            self.name = self.splited_data[4].replace('\n', '')
            self.base_color = False
        elif len(self.splited_data) == 2:
            self.gene = self.splited_data[0].replace(' ', '')
            self.color = self.splited_data[1].replace(' ', '').replace('\n', '')
            self.base_color = True
        else:
            raise ValueError(f'Unknown domain format: {input}')


class DomainDataset:
    def __init__(self, domain_content_list: List[str]):
        self.domains = self._extract_all_domain_content(domain_content_list)
        self._index = 0 
        self._size = len(self.domains)
    
    def _extract_all_domain_content(self, domain_content_list: List[str]) -> List['DomainEntity']:
        domains = []
        
        for content in domain_content_list.split('\n'):
            try:
                # Process each line in the content string
                for line in content.splitlines():
                    # Ignore comments and empty lines
                    if line and line[0] != '#' and line.strip():
                        domains.append(DomainEntity(line))
            except Exception as e:
                print(f'Domain_Dataset error: {e}')
        
        return domains

    def __len__(self):
        return len(self.domains)

    def __iter__(self) -> Iterator['DomainEntity']:
        self._index = 0  # Reset index for new iteration
        return self
    
    def __next__(self) -> 'DomainEntity':
        if self._index < self._size:
            domain = self.domains[self._index]
            self._index += 1
            return domain
        else:
            raise StopIteration

    def filter_by_fasta(self, FastaDataset: 'FastaDataset') -> 'DomainDataset':
        filtered_domains = []
        for domain in self.domains:
            for fasta in FastaDataset:
                if domain.gene == fasta.prot_gene:
                    filtered_domains.append(domain)
                    break

        self.domains = filtered_domains
        return self
