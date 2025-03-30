import colorsys
from dataclasses import dataclass
from typing import List, Tuple, Optional
import os
import sys
import copy
import re

from pycirclize import Circos as circos

# from matplotlib_venn import venn2
# from matplotlib_venn import venn3
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

from xldg.core import CrossLinkDataset, FastaDataset, DomainDataset


@dataclass
class CircosConfig:
    # File input
    fasta: FastaDataset
    domains: Optional[DomainDataset] = None
    
    # Text input
    legend: Optional[str] = None
    title: Optional[str] = None
    
    # Figure configs
    figsize: Tuple[float, float] = (9, 9)
    label_interval: int = 20
    space_between_sectors: int = 5
    domain_legend_distance: float = 1.25
    xl_legend_distance: float = 1.3
    xl_counter_distance: float = -0.15
    legend_distance: float = -0.15
    
    # Font configs
    title_font_size: int = 14
    ruler_font_size: int = 14
    legend_font_size: int = 14
    prot_font_size: int = 14
    
    # Figure elements plotting configs
    plot_all_proteins: bool = False
    plot_protein_ids: bool = True
    plot_counter: bool = True
    plot_xl_legend: bool = True
    plot_domain_legend: bool = True
    
    # XL configs
    min_rep: int = 1
    max_rep: int = sys.maxsize
    plot_interprotein_xls: bool = True
    plot_intraprotein_xls: bool = True
    plot_homotypical_xls: bool = True


class Circos:  
    def __init__(self, xls: 'CrossLinkDataset', config: 'CircosConfig'):
        self.config = copy.deepcopy(config)
        self.xls = copy.deepcopy(xls)

        self.xls.filter_by_replica(self.config.min_rep, self.config.max_rep)

        if self.config.plot_interprotein_xls is False:
            self.xls.remove_interprotein_crosslinks()
        if self.config.plot_intraprotein_xls is False:
            self.xls.remove_intraprotein_crosslinks()
        if self.config.plot_homotypical_xls is False:
            self.xls.remove_homotypic_crosslinks()


        self.fasta = copy.deepcopy(config.fasta)
        if config.plot_all_proteins is False:
            self.fasta.filter_by_crosslinks(self.xls)
        
        self.domains = None
        if self.config.domains is not None:
            self.domains = copy.deepcopy(self.config.domains)
            self.domains.filter_by_fasta(self.fasta)

        self.fig = None
        
        self.sectors = {prot.prot_gene: prot.seq_length for prot in self.fasta}
        self.prot_colors = self._assign_colors()
        self.circos = circos(self.sectors, space=self.config.space_between_sectors)

        # XL colors
        self.heterotypic_intraprotein_xl_color = '#21a2ed' # Blue
        self.heterotypic_interprotein_xl_color ='#00008B' # Dark Blue
        self.homotypic_xl_color = '#ed2b21' # Red
        self.general_xl_color = '#7d8082' # Grey
     
    def save(self, path: str) -> None:
        if len(self.xls) == 0:
            print(f'WARNING: No CrossLinkEntitys detected! Aborted save to {path}')
            return

        folder_path = os.path.dirname(path)
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
        
        self._plot_sectors()
        self._plot_xls()
        
        if (self.config.legend is not None):
            self._plot_user_legend()

        if (self.config.domains is not None and len(self.config.domains) != 0):
            self._plot_domains()
            
        if self.config.plot_xl_legend:
            self._plot_xl_legend()
        
        if self.config.plot_counter:
            self._plot_counter()
        
        if (self.config.title is not None):
            self._plot_title()

        self.fig.savefig(path)
        plt.close(self.fig)
        print(f'Circos plot saved to {path}')
    

    def _assign_colors(self) -> None:
        prot_colors = {}
        i = 0
        if self.domains is None:
            length = len(self.sectors)
            new_colors = self._generate_summer_colors(length)
            for prot in self.sectors:
                prot_colors[prot] = new_colors[i]
                i += 1
        else:
            for prot in self.sectors:
                prot_colors[prot] = '#C0C0C0'
                for domain in self.domains:
                    if domain.base_color is False:
                        continue
                    if prot == domain.gene:
                        prot_colors[prot] = domain.color
                        break        
                
        return prot_colors
    
    def _generate_summer_colors(self, num_colors: int) -> List[str]:
         summer_colors = []
         hue = 0.0  # Start at red (Hue 0)

         # Generate summer colors
         for _ in range(num_colors):
             lightness = 0.7  # High lightness for vibrant colors
             saturation = 0.8  # High saturation for bright colors
             r, g, b = colorsys.hls_to_rgb(hue, lightness, saturation)
             r = int(r * 255)
             g = int(g * 255)
             b = int(b * 255)
             summer_colors.append(f'#{r:02x}{g:02x}{b:02x}')

             # Increment hue for next color (golden ratio to avoid repeating colors)
             hue = (hue + 0.618033988749895) % 1.0

             # Skip hues to focus on typical summer colors (yellow, green, blue, pink)
             if 0.15 < hue < 0.3 or 0.55 < hue < 0.7:
                 hue = (hue + 0.2) % 1.0

         return summer_colors
    
    def _plot_sectors(self) -> None:
        for sector in self.circos.sectors:
            track = sector.add_track((92, 100))
            track.axis(fc = self.prot_colors[sector.name])
            if self.config.plot_protein_ids:
                sector.text(sector.name, color = '#3A3B3C', r = 110, size = self.config.prot_font_size)

            if self.domains != None:
                for domain in self.domains:
                    if domain.gene != sector._name or domain.base_color is True:
                        continue
                    track2 = sector.add_track((92, 100))
                    track2.rect(domain.start, domain.end, fc=domain.color)
            
            track._start += 1 # Remove zero lable of the plot
            track.xticks_by_interval(self.config.lable_interval, label_size = self.config.ruler_font_size)
            track._start -= 1

    def _plot_xls(self) -> None:
        for xl, site_count in self.xls.xls_site_count.items():
            xl_color = self.heterotypic_intraprotein_xl_color
            plane = 2

            protein_1 = self.fasta.find_gene_by_fasta_header(xl.protein_1)
            protein_2 = self.fasta.find_gene_by_fasta_header(xl.protein_2)
            if protein_1 == None or protein_2 == None:
                continue

            if xl.is_homotypical:
                xl_color = self.homotypic_xl_color
                plane = 3
            elif xl.is_interprotein:
                xl_color = self.heterotypic_interprotein_xl_color
            
            self.circos.link((protein_1, xl.num_site_1, xl.num_site_1), (protein_2, xl.num_site_2, xl.num_site_2), ec=xl_color, zorder=plane, lw=site_count)
        
        self.fig = self.circos.plotfig(figsize = self.config.figsize)
    
    def _plot_counter(self) -> None:
        total_xls_sites = 0
        site_counter = {}
        
        for xl, site_count in self.xls.xls_site_count.items():
            protein_1 = self.fasta.find_gene_by_fasta_header(xl.protein_1)
            protein_2 = self.fasta.find_gene_by_fasta_header(xl.protein_2)
            if protein_1 == None or protein_2 == None:
                continue
            
            total_xls_sites += 1
            
            if site_count in site_counter:
                site_counter[site_count] += 1
            else:
                site_counter[site_count] = 1
                
        sorted_site_counter = dict(sorted(site_counter.items()))
        
        if total_xls_sites > 0:
            text_lable = f'Total unique XLs: {total_xls_sites}\n'
            
            for key, value in sorted_site_counter.items():
                ending = ''
                if key > 1:
                    ending = 's'
                    
                text_lable += f'{key} replica{ending} unique XLs: {value}\n'
            
            self.fig.text(self.config.xl_counter_distance, 0.98, text_lable, fontsize=self.config.legend_font_size, va='top', ha='left')
      
    def _plot_user_legend(self) -> None:
        if self.config.legend != None:
            self.fig.text(self.config.legend_distance, 0.00, self.config.legend, va='bottom', ha='left', fontsize=self.config.legend_font_size)
           
    def _plot_domains(self) -> None:
        domains = [
            {'color': domain.color, 'label': domain.name}
            for domain in self.domains
            if domain.base_color is False
        ]
        legend_patches = []
        reference_buffer = []
        for item in domains:
            reference = item['color'] + item['label']
            if(reference in reference_buffer):
                continue
            
            check = item['label'].replace(' ', '')
            if(check != ''):
                legend_patches.append(mpatches.Patch(facecolor=item['color'], label=item['label'], linewidth=0.5, edgecolor='#3A3B3C'))
                reference_buffer.append(reference)
        
        if self.config.plot_domain_legend is True and len(legend_patches) != 0:
            self.fig.legend(handles=legend_patches, loc='lower right', bbox_to_anchor=(self.config.domain_legend_distance, 0), fontsize=self.config.legend_font_size)
    
    def _plot_xl_legend(self) -> None:
        most_frequent_xl = 0
        exhist_interprotein_xl = False
        exhist_intraprotein_xl = False
        exhist_homotypcal_xl = False

        for xl, site_count in self.xls.xls_site_count.items():
            if most_frequent_xl < site_count:
                most_frequent_xl = site_count

            if xl.is_homotypical:
                exhist_homotypcal_xl = True
            elif xl.is_interprotein:
                exhist_interprotein_xl = True
            else:
                exhist_intraprotein_xl = True
                
            
        if most_frequent_xl == 0:
            return

        legend_info = []
        if exhist_intraprotein_xl is True and self.config.plot_intraprotein_xls is True:
            legend_info.append({'label': 'Intraprotein unique XLs', 'color': self.heterotypic_intraprotein_xl_color, 'linewidth': 2})

        if exhist_interprotein_xl is True and self.config.plot_interprotein_xls is True:
            legend_info.append({'label': 'Interprotein unique XLs', 'color': self.heterotypic_interprotein_xl_color, 'linewidth': 2}) 

        if exhist_homotypcal_xl is True and self.config.plot_homotypical_xls is True:
            legend_info.append({'label': 'Homotypic unique XLs', 'color': self.homotypic_xl_color, 'linewidth': 2})

        if self.config.min_rep == 1:
            legend_info.append({'label': '1-replica unique XLs', 'color': self.general_xl_color, 'linewidth': 1})
        
        if most_frequent_xl > 1:
            for i in range(2, most_frequent_xl + 1):
                if i < self.config.min_rep:
                    continue

                legend_info.append({'label': f'{i}-replicas unique XLs', 'color': self.general_xl_color, 'linewidth': i}) 
        
        legend_handles = [Line2D([0], [0], color=info['color'], linewidth=info['linewidth'], label=info['label']) for info in legend_info]
        self.fig.legend(handles=legend_handles, loc='upper right', bbox_to_anchor=(self.config.xl_legend_distance, 1), fontsize=self.config.legend_font_size)
    
    def _plot_title(self) -> None:
        if self.config.title is not None:    
            self.fig.text(0.5, 1.05, self.config.title, ha='center', va='center', fontsize=self.config.title_font_size)
    
    def set_xls_colors(self, 
                      heterotypic_intraprotein_xl_color = '#21a2ed', 
                      heterotypic_interprotein_xl_color = '#00008B', 
                      homotypic_xl_color = '#ed2b21', 
                      general_xl_color = '#7d8082') -> None:

        def _valid_hex_color(color) -> str:
            hex_color_pattern = r'^#([0-9A-Fa-f]{3}){1,2}$'
            if re.match(hex_color_pattern, color):
                return color
            raise ValueError(f'ERROR! Invalid hex color: {color}')

        self.heterotypic_intraprotein_xl_color = _valid_hex_color(heterotypic_intraprotein_xl_color)
        self.heterotypic_interprotein_xl_color = _valid_hex_color(heterotypic_interprotein_xl_color)
        self.homotypic_xl_color = _valid_hex_color(homotypic_xl_color)
        self.general_xl_color = _valid_hex_color(general_xl_color)


# # @dataclass
# # class VennConfig:
# #     labels: List[str]
# #     domain_colors: List[str]
# #     overlap_colors: List[str]
# #     legend: Optional[str] = None
# #     title: Optional[str] = None
# #     title_font: int = 16
# #     legend_font: int = 16
# #     plot_legend: bool = True
# #     dpi: int = 600
# #     figsize: Tuple[float, float] = (9, 9)

# # # xls: 'CrossLinkDataset'

# # # TODO: test and improve
# # def build_ven2_diagram_of_xl_sites(save_path: str, xls_list1: 'CrossLinkDataset', xls_list2: 'CrossLinkDataset', config: 'Venn2_Config') -> None:
# #     self.config = copy.deepcopy(config)
# #     plt.figure(figsize=(10, 10), dpi=600)
    
# #     set1 = set([str(sublist) for sublist in xls_list1])
# #     set2 = set([str(sublist) for sublist in xls_list2])
# #     venn = venn2([set1, set2], (label1, label2))

# #     # Customize colors
# #     venn.get_patch_by_id('10').set_color('#9AE66E') # pastel green
# #     venn.get_patch_by_id('01').set_color('#FAF278') # pastel yellow
# #     venn.get_patch_by_id('11').set_color('#87D5F8') # pastel blue

# #     # Label the regions with the number of elements
# #     for subset in ('10', '01', '11'):
# #         if venn.get_label_by_id(subset):
# #             venn.get_label_by_id(subset).set_text(f'{venn.get_label_by_id(subset).get_text()}')

# #     # Customize font size
# #     for text in venn.set_labels:
# #         text.set_fontsize(size)

# #     for text in venn.subset_labels:
# #         if text:  # Check if the subset label is not None
# #             text.set_fontsize(size)
# #     if title != None:
# #         plt.title(title).set_fontsize(18)

# #     plt.savefig(save_path)
# #     print(f'Venn2 diagram saved to {save_path}')
