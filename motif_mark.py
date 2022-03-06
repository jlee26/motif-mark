#!/usr/bin/env python 

import Bioinfo
import argparse
import re
import cairo
import math
#from IPython.display import SVG, display, Image

# import sys

### Goal: With a list of motifs in a file and a fasta file with sequences, create a figure (1), 
###         and within the figure, show the motifs (all different kinds) by color coding.
###         Have a lened (which motif with which color). Title as the header? -- somehow indicate which graph goes with which sequence.
###         Show additional information, such as exons, start and end position?, 
### Pseudocode
### 1. Create OOP:
###     - For all the motifs
###     - For all the sequences (from FASTA)
###     - For all the introns (or exons)
###     - For all the genes?
### 2. Want to read the fasta file and make all the sequence in one line. Find the length of the sequence.
### 3. Find motifs in the sequence (regex). If found, keep the start & end position. Move to next motifs
### 4. Find exons in the sequence. If found, keep the start & end position.
### 5. Visualization:
### 6. Create a ____ dimension blank figure.
### 7. Createa sub-figure title at the top.
### 8. Create a horizontal line, length of the sequence. Find the start and end position of the motif and draw a colored box.
### 9. Find start and end position for exon and draw a black box.

### QUESTION: Can i require a GTF file so that we can say where the introns/exons are?
### QUESTION: Do we have to make the overlapped motifs staggered? Can we annotate the start position aabove?
### Use UNIX command line to only find exons?

# multiple fasta files?

#def get_args():
#    parser = argparse.ArgumentParser(description="Files")
#    parser.add_argument("-f", "--file", help="FASTA file", required=True)
#    parser.add_argument()
#    return parser.parse_args()
#args = get_args()

class Motif:
    '''This class represents all the motifs. It will take in a single text file that contains all the motifs in each line (not case-sensitive).'''
    def __init__(self, file):
        '''DOCUMENT'''
        ## Data ##
        self.file = file
    ## Methods ##
    def extract_motif(self):
        '''Creates a list with all the motifs in lower case. Since it's a list, may be good to not have too much motifs.'''
        with open(self.file, "r") as f:
            count = 0
            dict_motif = {}
            for line in f.readlines():
                dict_motif[count] = line.lower().strip()
                count += 1
                #line = f.read().lower().splitlines()
        return dict_motif
    def y_converter(self):
        '''Y can be cytosine or thymine (or uracil)'''
        dict_motifs = self.extract_motif()
        dict_possible_motifs = {}
        for key,motif in dict_motifs.items():
            if "y" in motif or "u" in motif:
                convert = motif.replace("y","[c or t]").replace("u","t")
                dict_possible_motifs[motif] = convert
                #print(motif)
            else:
                dict_possible_motifs[motif] = motif
        return dict_possible_motifs
         
#motif_findings = Motif('Fig_1_motifs.txt')
#motif_findings.extract_motif()
#motif_findings.y_converter()

class Exon:
    ''''''
    def __init__(self, file):
        '''DOCUMENT'''
        ## Data ##
        self.file = file
    
    ### Methods ###
    ### zcat Homo_sapiens.GRCh38.105.gtf.gz | awk '$3 ~ /exon/' > Homo_sapiens.GRCh38.105_exons.gtf
    def exon_pos(self):
        dict_exon_pos = {}
        with open(self.file, "r") as f:
            for line in f.readlines():
                line = line.strip()
                exon_id = re.search(r'(?<=exon_id ")[^"\s]*',line)
                line = line.split()
                chrom = line[0]
                start = line[3]
                end = line[4]
                if chrom not in dict_exon_pos.keys():
                    dict_exon_pos[chrom] = {}
                    dict_exon_pos[chrom][exon_id.group()] = start, end
                else:
                    dict_exon_pos[chrom][exon_id.group()] = start,end
        return dict_exon_pos
# test
# exon_findings = Exon("exons_sample.gtf")

class Sequence:
    '''This class represents all the sequences in a FASTA file.'''
    def __init__(self, file):
        '''DOCUMENT'''
        ## Data ##
        self.file = file
        self.motif = Motif('Fig_1_motifs.txt')
        self.exon = Exon("Homo_sapiens.GRCh38.105_exons.gtf")

    ### Methods ##
    def extract_sequence(self):
        '''Extract all the sequence from FASTA file.'''
        Bioinfo.oneline_fasta(self.file, "oneline_fasta.txt")
    def dictionary_sequence(self):
        ''''''
        self.extract_sequence()
        with open("oneline_fasta.txt", "r") as f:
            count = 0
            dict_seq = {}
            header= ""
            for line in f.readlines():
                if line[0] == ">":
                    header = line.strip()
                    #print(header)
                else:
                    dict_seq[header] = line.strip()
                    count += 1
        return dict_seq
    def seq_length(self):
        dict_seq = self.dictionary_sequence()
        dict_length = {}
        for header,seq in dict_seq.items():
            dict_length[header] = len(seq)
        return dict_length
    def find_motif(self):
        dict_possible_motifs = self.motif.y_converter()
        dict_seq = self.dictionary_sequence()
        dict_motif_in_seq = {header: {} for header,seq in dict_seq.items()}
            ##### create dictionary, named after header for every sequence (header). key = start pos, assuming there is not overlap between different motifs. value = motifs. Will be used for visualization
        for header,seq in dict_seq.items():
            for key,motifs in dict_possible_motifs.items():
                #print(header)
                for match in re.finditer(rf'{motifs}', seq):
                    dict_motif_in_seq[header][match.start()] = key
                    #print(key, match.start(), match.end())
        return dict_motif_in_seq
    def find_exon(self):
        dict_exon_pos = self.exon.exon_pos()
        dict_length = self.seq_length()
        for header,length in dict_length.items():
            chrom = re.search(r'(?<=chr")[^:\s]*',header)
            
        print(header)

        

# sequence_findings = Sequence('Figure_1.fasta')
# sequence_findings.extract_sequence()
# sequence_findings.dictionary_sequence()
# sequence_findings.find_motif()

class Visualization:
    ''''''
    def __init__(self):
        '''DOCUMENT'''
        ## Data ##
        self.sequence = Sequence('Figure_1.fasta')
        self.motif = Motif('Fig_1_motifs.txt')

    ### Methods ##
    #def show_img(self):
    #    ''''''
    #    display(Image(filename=file))
    #def show_svg(self):
    #    ''''''
    #    display(SVG(filename=file))
    def motif_color(self):
        dict_possible_motifs = self.motif.y_converter()
        dict_motif_color = {}
        list_color_temp = [[193, 187, 221], [255,189,230],[252, 252, 148],[166,205,229]] #purple, pink, yellow, blue
        list_color_temp_01 = []
        for col in list_color_temp:
            list_color_temp_01.append([x/255 for x in col])
        #print(list_color_temp_01)
        count = 0
        for motif,converter in dict_possible_motifs.items():
            dict_motif_color[motif] = list_color_temp_01[count]
            count += 1
        return dict_motif_color
    def figure(self):
        dict_seq = self.sequence.dictionary_sequence()
        dict_length = self.sequence.seq_length()
        dict_motif_in_seq = self.sequence.find_motif()
        dict_motif_color = self.motif_color()
        max_length = dict_length[max(dict_length, key=dict_length.get)]
        #print(max_length)
        width = max_length + 50 + 300 #50 is for buffer? 300 is for legend space (250 legend, 50 buffer)
        tot_read = len(self.sequence.dictionary_sequence())
        height = tot_read*300

        #print(max_length)
        surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width, height)
        context = cairo.Context(surface)

        #background color
        context.rectangle(0,0, width,height)
        context.set_source_rgb(1,1,1)
        context.fill()

        #legend box
        context.rectangle(max_length + 50, 75, 250, 250)
        context.fill_preserve()
        context.set_source_rgb(0,0,0)
        context.set_line_width(2)
        context.stroke()
        context.set_font_size(18)
        context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL,cairo.FONT_WEIGHT_NORMAL)
        context.move_to(max_length + 50 + 100, 75+25)
        context.show_text("Legend")
        

        #legend
        count=0
        for motif,color in dict_motif_color.items():
            context.set_source_rgb(0,0,0)
            context.set_font_size(16)
            context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL,cairo.FONT_WEIGHT_NORMAL)
            context.move_to(max_length+50+25+25,75+25+40+(count*50))
            context.show_text(motif)

            
            context.set_source_rgb(dict_motif_color.get(motif)[0],dict_motif_color.get(motif)[1],dict_motif_color.get(motif)[2])
            context.rectangle(max_length+50+25,75+25+30+(count*50),len(motif),15)
            context.fill()


            count +=1

        count=0
        for header,length in dict_length.items():
            #header
            context.set_source_rgb(0,0,0)
            context.set_font_size(20)
            context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL,cairo.FONT_WEIGHT_NORMAL)
            context.move_to(25+10,75+(count*300))
            context.show_text(header)

            #sequence (line)
            context.set_line_width(1)
            context.move_to(25, 150+(count*300))
            context.line_to(25+length,150+(count*300))
            context.set_source_rgb(0,0,0)
            context.stroke()

            #motifs (rectangle)
            ##### use the dictionary from Sequence. based on dictionary name, plot on sequence. if value is certain motif, use certain color and create rect.
            for h,pairs in dict_motif_in_seq.items():
                #print(h, pairs)

                if h == header:
                    for start_pos,motif in pairs.items():
                        context.set_source_rgb(dict_motif_color.get(motif)[0],dict_motif_color.get(motif)[1],dict_motif_color.get(motif)[2])
                        context.rectangle(25+start_pos,100+(count*300),len(motif),100)
                        context.fill()
                else:
                    continue
            
            #exons (rectangle)
            #context.rectangle(25+10,25-10,75-10-35,10+10)
            #context.fill()

            count += 1

        surface.write_to_png("Figure_motifs.png")


figure_visualization = Visualization()
figure_visualization.figure()
