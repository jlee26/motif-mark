#!/usr/bin/env python 

import Bioinfo
import argparse
import re
import cairo
import math
import random

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

def get_args():
   parser = argparse.ArgumentParser(description="Files")
   parser.add_argument("-f", "--fasta", help="Name of the FASTA file. Ex: 'Figure_1.fasta'", required=True)
   parser.add_argument("-m", "--motif", help="Name of the motif file. Ex: 'Fig_1_motifs.txt'", required=True)
   return parser.parse_args()
args = get_args()

class Motif:
    '''This class represents all the motifs. 
    It will take in a single text file that contains all the motifs on each line.
    These motifs is NOT case-sensitive (motifs can be both capitalized, lowercaesd, or mixture of both; these are treated the same).'''
    def __init__(self, file):
        '''file : this is the motif text file.'''
        ## Data ##
        self.file = file
    ## Methods ##
    def extract_motif(self):
        '''Creates a dictioanry with all the motifs in lower case. Key:count, value:motif(lowercase)'''
        with open(self.file, "r") as f:
            count = 0
            dict_motif = {}
            for line in f.readlines():
                dict_motif[count] = line.lower().strip()
                count += 1
                #line = f.read().lower().splitlines()
        return dict_motif
    def y_converter(self):
        '''Y's are treated as anonymous nucleotides and is converted to cytosine (C) or thymine (T).
        U's are converted to T's.
        Creates a dictionary with converted nucleotides. Key:motifs with anonymous nucleotides, value:motifs with converted nucleotides.'''
        dict_motifs = self.extract_motif()
        dict_possible_motifs = {}
        dict_anonymous = {"u":"t", "w":"[a or t]", "s":"[c or g]", "m":"[a or c]", "k":"[g or t]", "r":"[a or g]","y":"[c or t]", "n":"[a or c or g or t]"}
        for key,motif in dict_motifs.items():
            for symbol,rep in dict_anonymous.items():
                if symbol in motif and motif in dict_possible_motifs.keys():
                    
                    val = dict_possible_motifs[motif]
                    convert = val.replace(symbol,rep)
                    dict_possible_motifs[motif] = convert
                elif symbol in motif and motif not in dict_possible_motifs.keys():
                    convert = motif.replace(symbol,rep)
                    dict_possible_motifs[motif] = convert
                else:
                    if motif in dict_possible_motifs.keys():
                        continue
                    else:
                        dict_possible_motifs[motif] = motif
        return dict_possible_motifs
         
#motif_findings = Motif('Fig_1_motifs.txt')
#motif_findings.extract_motif()
#motif_findings.y_converter()

class Sequence:
    '''This class represents all the sequences in a FASTA file. This class will interact with Motif Class.
    This class can extract the sequences from the FASTA so that the multi-lined sequence is converted to one-line.
    This class can create a dictionary, where the values are the sequences. There is 2 options, all lower case or as is. Captializatiosn are the exons.
    This class can find the length of the sequnce.
    Thsi class can find the motifs start position. 2 options are available: unsorted and sorted.'''
    def __init__(self, file):
        '''file : this is a fasta file.'''
        ## Data ##
        self.file = file
        self.motif = Motif(args.motif)
        #self.exon = Exon("Homo_sapiens.GRCh38.105_exons.gtf")

    ### Methods ##
    def extract_sequence(self):
        '''Extract multi-lined sequence from FASTA file and put into one line.
        Creates a new file named 'oneline_fasta.txt' that has one line of sequence.'''
        Bioinfo.oneline_fasta(self.file, "oneline_fasta.txt")
    def dictionary_sequence(self):
        '''Returns a dictionary with DNA sequences, all lowercase.
        This is lower case to help find motifs using this sequence.
        Key: header of the FASTA file, Value: DNA sequence'''
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
                    dict_seq[header] = line.strip().lower()
                    count += 1
        return dict_seq
    def dictionary_sequence_upper(self):
        '''Returns a DNA sequence from FASTA file, as is with Upper cases.
        The upper cases indicates EXONS. Used to visualize exons/introns.
        Key: header of the FASTA file, Value: DNA sequence from FASTA (with upper cases)'''
        self.extract_sequence()
        with open("oneline_fasta.txt", "r") as f:
            count = 0
            dict_seq_upper = {}
            header= ""
            for line in f.readlines():
                if line[0] == ">":
                    header = line.strip()
                    #print(header)
                else:
                    dict_seq_upper[header] = line.strip()
                    count += 1
        return dict_seq_upper 
    def seq_length(self):
        '''Returns a dictionary with the length of each sequence.
        Can be used to determine longest length and set the visualization figure dimensions.
        Key: header of the FASTA file, Value: length of the sequence'''
        dict_seq = self.dictionary_sequence()
        dict_length = {}
        for header,seq in dict_seq.items():
            dict_length[header] = len(seq)
        return dict_length
    def find_motif(self):
        '''Returns a dictionary of the start location of the motifs within a dictionary of headers.
        Key: header of the FASTA file, Value:{Key: start position of the motifs, Value: motif sequence}'''
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
    def sort_dict_motif_in_seq(self):
        '''Returns a sorted dictioanry of function find_motif.
        Sorted dictionary of the start location of themotifs within a dictionary of headers.
        This is used for visualization of the motifs, so that the motifs are staggered correctly.'''
        dict_motif_in_seq = self.find_motif()
        dict_sorted_motif_in_seq = {header: {} for header,motif in dict_motif_in_seq.items()}
        for header,motif in dict_motif_in_seq.items():
            dict_sorted_motif_in_seq[header] = dict(sorted(dict_motif_in_seq[header].items()))
        return dict_sorted_motif_in_seq

# sequence_findings = Sequence('Figure_1.fasta')
# sequence_findings.extract_sequence()
# sequence_findings.dictionary_sequence()
# sequence_findings.find_motif()
# sequence_findings.sort_dict_motif_in_seq()


class Exon:
    '''This class contains all the exons. This class interacts with the Sequence Class.
    This can find the start and end position of the exons.'''
    def __init__(self):
        '''This class interacts with Sequence Class.'''
        ## Data ##
        self.sequence = Sequence(args.fasta)
    
    ### Methods ###
    ### zcat Homo_sapiens.GRCh38.105.gtf.gz | awk '$3 ~ /exon/' > Homo_sapiens.GRCh38.105_exons.gtf
    def exon_pos(self):
        '''This returns a dictioanry with start and end position of the motifs in the sequence.
        Key: header of the FASTA file, Value: [start position, end position of the motif]'''
        dict_seq_upper = self.sequence.dictionary_sequence_upper()
        #dict_exon_pos = {}
        dict_exon_pos = {header: {} for header,seq in dict_seq_upper.items()}

        #print(dict_seq)
        for header,seq in dict_seq_upper.items():
            count = 0
            
            for match in re.finditer(r'[A-Z]+',seq):
                #print(header, match.group(), match.start(), match.end())
                if header in dict_exon_pos.keys():
                    dict_exon_pos[header][count] = [match.start(), match.end()]
                    count += 1
                else:
                    continue
        #print(dict_exon_pos)
        return dict_exon_pos
           

# practice_exon = Exon()
# practice_exon.exon_pos()


class Visualization:
    '''This class creates visualization figure of the motifs. This interacts with the Sequence, Motif, and Exon class.
    This can find the color of the motifs.
    This can create the figure will all the reads with motifs, staggered.'''
    def __init__(self):
        '''This class interacts with sequence, motif, and exon Class.'''
        ## Data ##
        self.sequence = Sequence(args.fasta)
        self.motif = Motif(args.motif)
        self.exon = Exon()

    def motif_color(self):
        '''This function creates a dictioanry with all the motifs with color.
        Key: motifs, original without converting for anonymous nucleotides, Value: list of rgb'''
        dict_possible_motifs = self.motif.y_converter()
        dict_motif_color = {}
        #########this mannually creates color
        # list_color_temp = [[193, 187, 221], [255,189,230],[255, 215, 0],[166,205,229]] #purple, pink, yellow, blue
        # list_color_temp_01 = []
        # for col in list_color_temp:
        #     list_color_temp_01.append([x/255 for x in col])
        # #print(list_color_temp_01)
        # count = 0
        # for motif,converter in dict_possible_motifs.items():
        #     dict_motif_color[motif] = list_color_temp_01[count]
        #     count += 1
        ###########

        size = len(dict_possible_motifs)
        # for motif,converter in dict_possible_motifs.items():
        #     dict_motif_color[motif] = [random.sample(range(255),3)]
        list_color_temp = [random.sample(range(255),3) for i in range(size)]
        #print(list_color_temp)
        list_color_temp_01 = []
        for col in list_color_temp:
            list_color_temp_01.append([x/255 for x in col])
        #print(list_color_temp_01)
        count = 0
        for motif,converter in dict_possible_motifs.items():
            dict_motif_color[motif] = list_color_temp_01[count]
            count += 1
        #print(dict_motif_color)
        return dict_motif_color
    def figure(self):
        '''Creates a single figure. In each figure, there are many reads with motifs.'''
        dict_seq = self.sequence.dictionary_sequence()
        dict_length = self.sequence.seq_length()
        dict_sorted_motif_in_seq = self.sequence.sort_dict_motif_in_seq()
        dict_motif_color = self.motif_color()
        dict_exon_pos = self.exon.exon_pos()
        max_length = dict_length[max(dict_length, key=dict_length.get)]
        #print(max_length)
        width = max_length + 50 + 300 #50 is for buffer? 300 is for legend space (250 legend, 50 buffer)
        tot_read = len(self.sequence.dictionary_sequence())
        height = tot_read*300 + 75

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

        # title
        context.set_source_rgb(0,0,0)
        context.set_line_width(2)
        context.stroke()
        context.set_font_size(25)
        context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL,cairo.FONT_WEIGHT_NORMAL)
        context.move_to(max_length/2, 75)
        context.show_text("Motif Binding Sites")

        count=0
        overlap = 0
        for header,length in dict_length.items():
            #header
            context.set_source_rgb(0,0,0)
            context.set_font_size(20)
            context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL,cairo.FONT_WEIGHT_NORMAL)
            context.move_to(25+10,75+75+(count*250))
            context.show_text(header)

            #sequence (line)
            context.set_line_width(1)
            context.move_to(25, 150+(count*250)+75)
            context.line_to(25+length,150+(count*250)+75)
            context.set_source_rgb(0,0,0)
            context.stroke()

            #exons (rectangle)
            for h,exon_pos in dict_exon_pos.items():
                if h == header:
                    for i,pos in exon_pos.items():
                        start = pos[0]
                        end = pos[1]-pos[0]
                        context.set_source_rgb(0,0,0)
                        context.rectangle(25+start,125+(count*250)+75,end,50)
                        context.fill()
                else:
                    continue
            
            #motifs (rectangle)
            ##### use the dictionary from Sequence. based on dictionary name, plot on sequence. if value is certain motif, use certain color and create rect.
            for h,pairs in dict_sorted_motif_in_seq.items():
                
                #print(h, pairs)
                overlap = 0
                if h == header:
                    for start_pos,motif in pairs.items():
                        #print(start_pos)
                            
                            #overlap += 1
                        if overlap % 2 == 0:
                            context.set_source_rgb(dict_motif_color.get(motif)[0],dict_motif_color.get(motif)[1],dict_motif_color.get(motif)[2])
                            
                            context.rectangle(25+start_pos,120+(count*250)+75,len(motif),25)
                            context.fill()
                            overlap += 1
                            
                        elif overlap % 2 == 1:
                            context.set_source_rgb(dict_motif_color.get(motif)[0],dict_motif_color.get(motif)[1],dict_motif_color.get(motif)[2])

                            context.rectangle(25+start_pos,95+(count*250)+75,len(motif),25)
                            context.fill()
                            overlap += 1

                else:
                    continue
                


            count += 1
        
        #Figure info
        context.set_source_rgb(0,0,0)
        context.set_line_width(2)
        context.stroke()
        context.set_font_size(18)
        context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL,cairo.FONT_WEIGHT_NORMAL)
        context.move_to(50, height-100)
        context.show_text("Figure 1. The motif binding sites are shown for each read in a single FASTA file. The motifs are color coded")
        context.move_to(50, height-75)
        context.show_text("and are scaled to the sequence length. The motifs are staggered to ensure the motifs do not overlap.")
        context.move_to(50, height-50)
        context.show_text("The header of each visualizations indicate which chromosome and the start and end position of the sequence reads. ")
        context.move_to(50, height-25)
        context.show_text("Some visualizations will have a black rectangle \n in the center which indicates the exons and the horizontal line indicates introns.")


        size=len(args.fasta)
        filename = f"{args.fasta}"[:-6]
        surface.write_to_png(f"{filename}.png")


figure_visualization = Visualization()
figure_visualization.figure()
