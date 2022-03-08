#!/usr/bin/env python
#populate this file with any fuctions you have written so far
#How to import into working file: import Bioinfo

#CONSTANTS
DNA_BASES = set('ATGCatcg')
RNA_BASES = set('AUGCaucg')

#FUNCTIONS
def convert_phred(letter):
    """Converts a phred quality score from their single character into an integer"""
    return (ord(letter) -33)

if __name__ == "__main__":
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("convert_phred function is working!")

def qual_score(phred_score):
    '''Calculates the average quality score of the whole phred string (in characters).'''
    tot = 0
    for i in phred_score:
        iscore = convert_phred(i)
        tot += iscore
    size = len(phred_score)
    ave = tot/size
    return ave

if __name__ == "__main__":
    assert qual_score("FFHHHHHJJJJIJIJJJIJJJJJJIIIJJJEHJJJJJJJIJIDGEHIJJFIGGGHFGHGFFF@EEDE@C??DDDDDDD@CDDDDBBDDDBDBDD@") == 37.62105263157895, "wrong average for phred score"
    print("qual_score calculated the correct average quality score!")


def validate_base_seq(seq,RNAflag=False):
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    DNAbases = set('ATGCatcg')
    RNAbases = set('AUGCaucg')
    return set(seq)<=(RNAbases if RNAflag else DNAbases)

if __name__ == "__main__":
    assert validate_base_seq("AATAGAT") == True, "Validate base seq does not work on DNA"
    assert validate_base_seq("AAUAGAU", True) == True, "Validate base seq does not work on RNA"
    assert validate_base_seq("Hi there!") == False, "Validate base seq fails to recognize nonDNA"
    assert validate_base_seq("Hi there!", True) == False, "Validate base seq fails to recognize nonDNA"
    print("validate_base_seq works!")

def gc_content(DNA):
    '''Returns GC content of a DNA sequence as a decimal between 0 and 1.'''
    assert validate_base_seq(DNA), "String contains invalid characters"
    
    DNA = DNA.upper()
    Gs = DNA.count("G")
    Cs = DNA.count("C")
    return (Gs+Cs)/len(DNA)

if __name__ == "__main__":
    assert gc_content("GCGCGC") == 1
    assert gc_content("AATTATA") == 0
    assert gc_content("GCATGCAT") == 0.5
    print("gc_content correctly calculated GC content!")


def oneline_fasta(filename, outputfile):
    '''Takes two arguments: filename is the fasta file with multiple sequence lines and outputfile is the name of the file that will have a single line for the sequence.Returns an outputfile with the header and one line of sequence. Warning: make sure there are no pre-exisitng output file name, as this will overwrite any preexisitng file. Warning: the first line of the output file will be a blank line.'''
    with open(f"./{outputfile}", "w") as wf:
        with open(filename, "r") as rf:
            start = 0 
            for line in rf:
                if line.startswith(">") and start == 0:
                    wf.write(line)
                    start = 1
                elif line.startswith(">") and start == 1:
                    wf.write("\n"+line)
                else:
                    wf.write(line.strip())
