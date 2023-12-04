# from Bio.Blast import NCBIWWW
# from Bio.Blast import NCBIXML
import random
# # Define the DNA sequence
# dna_sequence = "TGGGCCTCATATTTATCCTATATACCATGTTCGTATGGTGGCGCGATGTTCTACGTGAATCCACGTTCGAAGGACATCATACCAAAGTCGTACAATTAGGACCTCGATATGGTTTTATTCTGTTTATCGTATCGGAGGTTATGTTCTTTTTTGCTCTTTTTCGGGCTTCTTCTCATTCTTCTTTGGCACCTACGGTAGAG"

# # Perform the BLAST search
# result_handle = NCBIWWW.qblast("blastn", "nr", dna_sequence)

# # Parse the BLAST results
# blast_records = NCBIXML.parse(result_handle)

# # Iterate over the blast records and print the species of each sequence
# for blast_record in blast_records:
#     for alignment in blast_record.alignments:
#         print("Title:", alignment.title)
#         print("Length:", alignment.length)
#         print("E-value:", alignment.hsps[0].expect)
#         print("Bit Score:", alignment.hsps[0].bit_score)
#         print("")
# try: 
#    f=open("dna2.fasta", 'r')
# except IOError:
#     print("file does not exist")

# seqs={}

# for line in f:
#     line=line.rstrip()
#     if line[0]=='>':
#         words=line.split()
#         name=words[0][1:]
#         seqs[name]=''
#     else:
#         seqs[name]=seqs[name]+line
      
           
     
# max_length = min(len(seq) for seq in seqs.values())
# count = sum(len(seq) == max_length for seq in seqs.values())

# print("Longest Sequence Length:", max_length)
# print("Count of Sequences with Longest Length:", count)
        
        
# # def get_extension1(filename):
# #     return(filename.split(".")[-1])

# # def get_extension2(filename):
# #     import os.path
# #     return(os.path.splitext(filename)[1])

# # def get_extension3(filename):
# #     return filename[filename.rfind('.'):][1:]
    

# # print(get_extension1('myfile.tar.gaz'))
# # print(get_extension2('myfile.tar.gaz'))
# # print(get_extension3('myfile'))



# from Bio.Seq import Seq

# def find_longest_orf_in_reading_frame_3(sequence):
#   """Finds the longest open reading frame (ORF) in a DNA sequence in reading frame 3.

#   Args:
#     sequence: A DNA sequence.

#   Returns:
#     A tuple containing the length of the longest ORF and the starting position of the longest ORF.
#   """

#   longest_orf_length = 0
#   longest_orf_start = 0
#   current_orf_length = 0
#   current_orf_start = 0
#   for i in range(2, len(sequence), 3):
#     codon = sequence[i:i + 3]
#     if codon == "ATG":
#       current_orf_start = i
#       current_orf_length = 3
#     elif codon in ["TAA", "TAG", "TGA"]:
#       if current_orf_length > longest_orf_length:
#         longest_orf_length = current_orf_length
#         longest_orf_start = current_orf_start
#       current_orf_length = 0
#       current_orf_start = 0
#     else:
#       current_orf_length += 3
#   return longest_orf_length, longest_orf_start

# def find_longest_orf_in_reading_frame_3_in_fasta_file(fasta_file):
#   """Finds the longest open reading frame (ORF) in reading frame 3 in any of the sequences in a FASTA file.

#   Args:
#     fasta_file: A FASTA file.

#   Returns:
#     A tuple containing the length of the longest ORF and the starting position of the longest ORF.
#   """

#   longest_orf_length = 0
#   longest_orf_start = 0
#   for sequence in parse_fasta_file(fasta_file):
#     orf_length, orf_start = find_longest_orf_in_reading_frame_3(sequence)
#     if orf_length > longest_orf_length:
#       longest_orf_length = orf_length
#       longest_orf_start = orf_start

#   return longest_orf_length, longest_orf_start

# def parse_fasta_file(fasta_file):
#   """Parses a FASTA file into a list of sequences.

#   Args:
#     fasta_file: A FASTA file.

#   Returns:
#     A list of sequences, where each sequence is a string.
#   """

#   sequences = []
#   current_sequence = ""
#   with open(fasta_file, "r") as f:
#     for line in f:
#       if line.startswith(">"):
#         if current_sequence:
#           sequences.append(current_sequence)
#         current_sequence = ""
#       else:
#         current_sequence += line.strip()
#   if current_sequence:
#     sequences.append(current_sequence)
#   return sequences

# if __name__ == "__main__":
#   fasta_file = "dna2.fasta"
#   longest_orf_length, longest_orf_start = find_longest_orf_in_reading_frame_3_in_fasta_file(fasta_file)

#   print("Longest ORF length:", longest_orf_length)


#   print("Starting position of the longest ORF:", longest_orf_start)



# from Bio.Seq import Seq


# def find_longest_orf_in_any_forward_reading_frame(sequence):
#   """Finds the longest open reading frame (ORF) in a DNA sequence in any forward reading frame.

#   Args:
#     sequence: A DNA sequence.

#   Returns:
#     The length of the longest ORF.
#   """

#   longest_orf_length = 0
#   for i in range(0, 3):
#     translated_sequence = sequence[i::3]
#     orf_length = len(max(translated_sequence.split("*"), key=len))
#     if orf_length > longest_orf_length:
#       longest_orf_length = orf_length

#   return longest_orf_length

# def find_longest_orf_in_any_forward_reading_frame_in_fasta_file(fasta_file):
#   """Finds the length of the longest open reading frame (ORF) appearing in any sequence and in any forward reading frame.

#   Args:
#     fasta_file: A FASTA file.

#   Returns:
#     The length of the longest ORF.
#   """

#   longest_orf_length = 0
#   for sequence in parse_fasta_file(fasta_file):
#     orf_length = find_longest_orf_in_any_forward_reading_frame(sequence)
#     if orf_length > longest_orf_length:
#       longest_orf_length = orf_length

#   return longest_orf_length

# def parse_fasta_file(fasta_file):
#   """Parses a FASTA file into a list of sequences.

#   Args:
#     fasta_file: A FASTA file.

#   Returns:
#     A list of sequences, where each sequence is a string.
#   """

#   sequences = []
#   current_sequence = ""
#   with open(fasta_file, "r") as f:
#     for line in f:
#       if line.startswith(">"):
#         if current_sequence:
#           sequences.append(current_sequence)
#         current_sequence = ""
#       else:
#         current_sequence += line.strip()
#   if current_sequence:
#     sequences.append(current_sequence)
#   return sequences

# if __name__ == "__main__":
#   fasta_file = "dna2.fasta"
#   longest_orf_length = find_longest_orf_in_any_forward_reading_frame_in_fasta_file(fasta_file)

#   print("Longest ORF length:", longest_orf_length)
  
  
  
  

# from Bio.Seq import Seq

# sequence = Seq("""GTCGATCGACACGACGCTCGCGCAGCGCGACGCGAAGGCCGCGTGAGCGCACGACGCGCGTCACACACCA
# CGAGCACAACGAACACGACCCCACTCTCACGGAGCCGACCATGGCCGACCTTCGCTGCACCATCGCGGGC
# ATCACTTCGCCGAACCCTTTCTGGCTGGCGTCCGCGCCGCCGACCGACAAGGCCTACAACGTGAACCGCG
# CGTTCGAGGCGGGCTGGGGCGGGGTCGTCTGGAAGACGCTCGGGCTCGATCCGCATGTCGTCAACGTCAG
# TTCGCGCTATGGCGCGGTGCAGTGGAACGGCCAGCGCATCGCGGGGCTGAACAACATCGAGCTGATCACC
# GACCGTCCGCTCGACGTGAACCTGAGAGAGATCGCGCAGGTGAAGCGCGACTGGCCGGACCGCGCGCTGA
# TCGTGTCGCTGATGGTGCCGTGCAACGAGCGCGACTGGAAATGGATCCTGCCGCTCGTCGAGGATACGGG
# CGCCGACGCGGTCGAGCTGAACTTCGGTTGTCCGCACGGGATGAGCGAGCGCGGGATGGGCGCGGCGGTC
# GGGCAGGTGCCCGAATATGTGGAGATGGTCACGCGCTGGGTGAAGGAAGGCACGAAGCTGCCGTGCCTCG
# TGAAGCTCACGCCGAACATCAGCGACATCCGGATGGGGTCGCGCGCCGCGTACAAGGGCGGCGCGGACGG
# CGTGTCGCTGATCAACACGATCAACTCGATCGTCGCGGTCGATCTCGACCATATGGCGCCGATGCCGACG
# GTCGACGGCAAGGGCACGCACGGCGGCTATTGCGGCCCGGCGGTCAAGCCGATCGCATTGAACATGGTCG
# CGGAGATCGCACGTGACCCGGAAACGCCGAACCTGCCGATCTCGGGCATCGGCGGCATCTCGTCATGGCG
# CGACGCGGCGGAGTTCATGGTGCTCGGCGCCGGCAGCGTGCAGGTGTGCACCGCCGCGATGCATTACGGA
# TTCCGGATCGTGTCGGACCTGGCCGACGGATTGTCGAACTGGATGGACGAGAAGGGCTACGCGACGCTCG
# ACGACATTCGCGGCCGCGCGGTGCCGAACGTGACCGACTGGAAATACCTGAACCTGAAATACGACATCAA
# GGCGCGTATCGACCAGGACCGCTGCATCCAGTGCGGGTTGTGCCATATCGCGTGCGAGGACACGTCGCAC
# CAGGCGATCACCGCGACGAAGGACGGCGTGCGGCATTTCGAAGTGGTCGATTCGGCGTGCGTCGGGTGCA
# ATCTTTGCATGCATGTGTGTCCGGTCGAGCAATGCATCACGATGGAGCGTGTCGATTCGGGCGACTACGC
# GAACTGGACCACGCATCCGAACAATCCGGCGAGCGCGGAGGCGGGGGCGAGTGCAGGCGCGGCGGCACCC
# GAGAAGCACGCGAAGAAGGCTGCTTGACGGCGTCCGGCGATGCGGGCCATCCTGCATCGCCGCCTTTCGT
# TCCACCCGGGCCGGCATCGAGTGATGCCGGCGTTGACGTTTTCGTGGAGTGAGTCAGATGAATCACGCAG
# CGAATCCCGCCGATCCCGATCGCGCCGCGGCGCAGGGCGGCAGCCTGTACAACGACGATCTCGCGCCGAC
# GACGCCGGCGCAGCGCACGTGGAAGTGGTATCACTTCGCGGCGCTGTGGGTCGGGATGGTGATGAACATC
# GCGTCGTACATGCTCGCGGCCGGGCTGATCCAGGAAGGCATGTCGCCGTGGCAGGCGGTGACGACGGTGC
# TGCTCGGCAACCTGATCGTGCTCGTGCCGATGCTGCTGATCGGCCATGCGGGCGCGAAGCACGGGATTCC
# GTACGCGGTGCTCGTGCGCGCGTCGTTCGGCACGCAGGGGGCGAAGCTGCCGGCGCTGCTGCGCGCGATC
# GTCGCGTGCGGCTGGTACGGGATCCAGACCTGGCTCGGCGGCAGCGCGATCTATACGCTGCTGAACATCC
# TGACCGGCAACGCGCTGCATGGCGCCGCGCTGCCGGTCATCGGCATCGGGTTCGGGCAGCTCGCATGCTT
# CCTCGTGTTCTGGGCGCTGCAGCTCTACTTCATCTGGCATGGCACCGATTCGATCCGCTGGCTCGAAAGC
# TGGTCGGCGCCGATCAAGGTCGTGATGTGCGTGGCGCTGGTGTGGTGGGCAACGTCGAAGGCGGGCGGCT
# TCGGCACGATGCTGTCGGCGCCGTCGCAGTTTGCCGCAGGCGGCAAGAAAGCCGGGCTGTTCTGGGCGAC
# CTTCTGGCCGGGGCTGACCGCGATGGTCGGCTTCTGGGCGACGCTCGCGCTGAACATCCCCGACTTCACG
# CGCTTCGCGCATTCGCAGCGCGACCAGGTGATCGGCCAGTCGATCGGGCTGCCGTTGCCGATGGCGCTGC
# TGTCGGTGGTGTCGGTCGTCGTGACGTCGGCGACCGTCGTGATCTACGGCAACGCGATCTGGGATCCGAT
# CGACCTGACGAGCCGGATGACGGGCATCGGCGTGGGCATCGCGCTCGTGATCCTCACGCTCGACACGATG
# TGCTGCAACCTCGCCGCGAATCTCGTCGGCCCGGCGTACGACTTCTCGAGCCTGTGGCCGAAGGCGATCT
# CGTACCGCACCGGCGGGATGATCACCGCGACGCTCGCGATCGTGATGATGCCGTGGAAGATCCTCGCGAC
# GACGGACGGCTACATCTTCACCTGGCTCGTCGGCTACTCGGCGCTGCTCGGGCCCGTGGCGGGGATCCTG
# ATGGTCGACTACTTCCTGATTCGCGGCACGCGGCTCGACACGCGCGCGCTGTTCGACGAGCGCGGCGGCT
# TCAGCTACGCGCGCGGCTGGAACCCGGCCGCGCTGGCCGCGCTCGCGGTCGGCGTGCTGCCGAACCTGCC
# CGGCTTCCTGCACACGGCGTTTCCGGCGTCGTTTCCGAACGTGCCGGCGTTCTTCAACACGCTTTACACG
# TACGCGTGGTTCGTCGGCCTCGTGCTGGCGTCATGCGTGTACGGCACCTGGATGAAGTGGCGCGCCGGAC
# AGCACGCGCAGATCGCGAGCGCCTGATTCGGCACCCGACAGTCAACGAGGAGGCAACCCCATGGCAATCC
# TGATTCGTGGCGGCACCGTGGTCGATGCGGACCGTTCCTACCGCGCGGACGTGCTCTGCGCAGCCCCGGA
# GGACGGCGGCACGATCCTGCAGATCGCCGGGCAGATCGATGCGCCGGCCGGCGCGACCGTCGTCGATGCG
# CACGACCAGTACGTGATGCCGGGCGGCATCGATCCGCATACGCACATGGAACTGCCGTTCATGGGCACGA
# CCGCGAGCGACGATTTCTACTCGGGTACGGCCGCCGGGCTCGCGGGCGGCACGACGAGCATCATCGACTT
# CGTGATCCCGAGCCCGAAGCAGCCGCTGATGGACGCGTTCCATGCCTGGCGCGGCTGGGCCGAGAAGGCG
# GCGGCCGACTACGGCTTCCACGTGGCCGTGACGTGGTGGGACGAGAGTGTGCACCGCGACATGGGCACGC
# TCGTGCGCGAACACGGCGTGTCGAGCTTCAAGCACTTCATGGCGTACAAGAACGCGATCATGGCCGACGA
# CGAGGTGCTCGTGAACAGCTTCTCGCGTTCGCTCGAACTCGGCGCGTTGCCGACCGTGCATGCGGAGAAC
# GGCGAGCTCGTGTTCCAGTTGCAGAAGGCGCTGCTCGCGCGCGGGATGACGGGGCCGGAGGCGCATCCGC
# TGTCGCGGCCGCCGGAGGTCGAGGGTGAGGCGGCGAATCGTGCGATCCGCATTGCGCAGGTGCTCGGCGT
# GCCGGTGTATATCGTGCATGTGTCCGCGAAGGACGCGGTCGATGCGATCACGAAGGCGCGCAGCGAAGGG
# CTGCGCGTGTTCGGCGAGGTGCTGCCGGGCCATCTGGTGATCGACGAGGCCGTCTATCGCGATCCGGACT
# GGACACGTGCGGCCGCGCACGTGATGAGCCCGCCGTTCCGCTCGGCCGAGCACCGCGAGGCGCTGTGGCG
# CGGGCTGCAGGCAGGGCAGCTGCATACGACGGCAACCGACCACTGCGTGTTCTGCGCGTCGCAGAAGGCG
# ATGGGCCGCGAGGATTTCACGAAGATCCCGAACGGCTGCGGCGGTGTCGAGGATCGCATGTCGGTGCTGT
# GGCATCACGGCGTGAATCATGGCCGCATCACGCCGAACGAGTTCGTGCGGATCACGTCGACGAACGCCGC
# GCAGATCTTCAACCTGTATCCGCGCAAGGGCGCCGTGCAGGTGGGCGCCGATGCCGACCTCGTCGTGTGG
# GACCCGGCCGCGACCAGGACGATCTCGGTGAAGACGCATCACCAGCAGGTCGATTTCAACGTGTTCGAGG
# GGATGACCGTACAAGGCGTCGCAACCCACACGCTCACGCGCGGCGCGCTCGCGTGGGCCGACGGCGATCT
# GCGTGCCGTGCGCGGCGCGGGCCGCTATCTGAAGCGCCCGCCGGCAGCCAGCTACTACGAGGCCGCGCGG
# ATCGCGAACCGGCTGCGCGAACCGCATCCGGTCGAGCGCGCCGGTTGAGCGTTGCGTATCGCGCGGGGCG
# TGTCGGTTCGAACGACACGCCCCGCGCATGTTTGAGCGTGCGTTTACGTTCGTGCCGGCACCGCTCGTGC
# CGCCTCTCCGCATCACGCCCATCCTCTCAATATTTGGGATGAATTGAGCGCGATCGCGCCTTGCCGATCT
# CCGGATACATAGAACAACTGAGCAAGTCGATGAAACACGCGATGTCGCGCAAATGCGACCATTTTGTTTG
# CGTTGTCGACGTGCATGCCGGCGAGTAATATCCACCGACGGCGT""")

# longest_orf_length = 0
# for i in range(0, 3):
#     orf_length = len(max(sequence[i::3].split("*"), key=len))
#     if orf_length > longest_orf_length:
#         longest_orf_length = orf_length

# print("Longest ORF length:", longest_orf_length)




# def find_most_frequent_repeat(sequences, repeat_length):
#   """Finds the most frequently occurring repeat of a given length in all sequences.

#   Args:
#     sequences: A list of DNA sequences.
#     repeat_length: The length of the repeats to find.

#   Returns:
#     A tuple containing the most frequent repeat and the number of times it occurs.
#   """

#   counts = {}
#   for sequence in sequences:
#     for i in range(len(sequence) - repeat_length + 1):
#       repeat = sequence[i:i + repeat_length]
#       if repeat not in counts:
#         counts[repeat] = 0
#       counts[repeat] += 1

#   most_frequent_repeat = max(counts, key=counts.get)
#   most_frequent_repeat_count = counts[most_frequent_repeat]

#   return most_frequent_repeat, most_frequent_repeat_count


# print(find_most_frequent_repeat(seqs.values(), 7))




# def find_longest_orf_in_reading_frame_2(sequence):
#   """Finds the longest open reading frame (ORF) in a DNA sequence in reading frame 2.

#   Args:
#     sequence: A DNA sequence.

#   Returns:
#     The length of the longest ORF.
#   """

#   longest_orf_length = 0
#   current_orf_length = 0
#   current_orf_start = 0
#   for i in range(1, len(sequence), 3):
#     codon = sequence[i:i + 3]
#     if codon == "ATG":
#       current_orf_start = i
#       current_orf_length = 3
#     elif codon in ["TAA", "TAG", "TGA"]:
#       if current_orf_length > longest_orf_length:
#         longest_orf_length = current_orf_length
#       current_orf_length = 0
#       current_orf_start = 0
#     else:
#       current_orf_length += 3
#   return longest_orf_length

# def find_length_of_longest_orf_appearing_in_reading_frame_2(sequences):
#   """Finds the length of the longest open reading frame (ORF) appearing in reading frame 2.

#   Args:
#     sequences: A list of DNA sequences.

#   Returns:
#     The length of the longest ORF appearing in reading frame 2.
#   """

#   longest_orf_length = 0
#   for sequence in sequences:
#     orf_length = find_longest_orf_in_reading_frame_2(sequence)
#     if orf_length > longest_orf_length:
#       longest_orf_length = orf_length

#   return longest_orf_length


# sequences = max_length = (seq for seq in seqs.values())

# length_of_longest_orf_appearing_in_reading_frame_2 = find_length_of_longest_orf_appearing_in_reading_frame_2(sequences)

# print("Length of the longest ORF appearing in reading frame 2:", length_of_longest_orf_appearing_in_reading_frame_2)




# from Bio import SeqIO
# record = SeqIO.read("dna2.fasta", 'fasta')
# table = 11
# min_pro_len = 100

# for strand, nuc in [(+1, record.seq), (-1, record.seq.reverse_complement())]:
#     for frame in range(3):
#         length = 3 * ((len(record)-frame) // 3) #Multiple of three
#         for pro in nuc[frame:frame+length].translate(table).split("*"):
#             if len(pro) >= min_pro_len:
#                 print("%s...%s - length %i, strand %i, frame %i" \
#                       % (pro[:30], pro[-3:], len(pro), strand, frame))
                
                
                
                
                
                
seq=''.join([random.choice('ATGC') for _ in range(5)])
print(seq)


def sequecne(s1,s2):
  i=0
  while i<len(s1) and i<len(s2) and s1[i]==s2[i]:
    i+=1
  return s1[:i]
print(sequecne('ATTGTGTATG','ATTGTGATACACAA'))



def complimnet(s):
  complementary={'A':'T', 'C':'G', 'G':'C', 'T':'A'}
  t=''
  for base in s:
    t=complementary[base]+t
  return t
print(complimnet('ATGC'))


# import requests

# url = "http://d28rh4a8wq0iu5.cloudfront.net/ads1/data/ERR266411_1.first1000.fastq"
# output_file = "file.fastq"

# response = requests.get(url)
# with open(output_file, "wb") as file:
#     file.write(response.content)

# print("File downloaded successfully.")

def readfast(file):
  sequences=[]
  qualities=[]
  
  with open(file) as fh:
    while True:
      fh.readline()
      seq =fh.readline().rstrip()
      fh.readline()
      qual=fh.readline().rstrip()
      if len(seq)==0:
        break
      sequences.append(seq)
      qualities.append(qual)
  return sequences,qualities

sequences,qualities=readfast("lambda_virus.fa")
print(qualities[:5])


# def pred33(qual):
#   return ord(qual)-33


# # print(pred33('#'))

  
      
# def createhist(qualities):
#   hist=[0]*50
#   for qual in qualities:
#     for phred in qual:
#       q=pred33(phred)
#       hist[q]+=1
#   return hist
# h=createhist(qualities)
# print(h)


# import matplotlib.pyplot as plt

  
# plt.bar(range(len(h)),h)
# plt.show()
      
    
# def findgc(reads):
#   gc=[0]*100
#   totals=[0]*100
#   for read in reads:
#     for i in range(len(read)):
#       if read[i]=='G' or read[i]=='C':
#         gc[i]+=1
#       totals[i]+=1
#   for i in range(len(gc)):
#     if totals[i]!=0:
#       gc[i]=gc[i]/float(totals[i])
#   return gc
# gc=findgc(qualities)
# plt.plot(range(len(gc)),gc)
# plt.show()


# import collections 
# count=collections.Counter()
# for seq in sequences:
#   count.update(seq)
# print(count)
  
  
    
def naive(p,t):
  occurnece=[]
  for i in range(len(t)-len(p)+1):
    match=True
    for j in range(len(p)):
      if t[j+i] != p[j]:
        match=False
        break
    if match:
        occurnece.append(i)
  return occurnece
t = "ABABDABACDABABCABAB"
p = "ABABC"
occurence=naive(p,t)
      
print(occurence)
        
t = ["ABABDABACDABABCABAB",'ads']
p = "ABABC"

import random
def generateRead(gene,numreads,readlen):
  reads=[]
  for _ in range(numreads):
    start=random.randint(0,len(gene)-readlen)-1
    reads.append(gene[start:start+readlen])
  return reads
reads=generateRead(sequences,100,100)



nummatched=0
for r in reads:
  matches=naive(r,sequences)
  if len(matches)>0:
    nummatched+=1
print('%d/%d reads matches exactly' %(nummatched,len(reads)))




def longest_prefix_suffix(string):
  """Returns the longest prefix of the given string that is also a suffix."""
  prefix_length = 0
  for i in range(len(string)):
    if string[i] == string[-i - 1]:
      prefix_length += 1
    else:
      break
  return prefix_length

string = "CACACTGCACAC"
prefix_length = longest_prefix_suffix(string)
print(prefix_length)





# import re

# def find_leftmost_occurrence_in_fasta(fasta_file, pattern):
#   """Returns the offset of the leftmost occurrence of the pattern in the FASTA file."""

#   leftmost_occurrence = -1
#   with open(fasta_file, 'r') as f:
#     for line in f:
#       if line.startswith('>'):
#         continue
#       else:
#         match = re.search(pattern, line)
#         if match:
#           leftmost_occurrence = match.start()
#           break
#   return leftmost_occurrence

# def reverse_complement(sequence):
#   """Returns the reverse complement of the sequence."""

#   complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
#   reversed_sequence = sequence[::-1]
#   reverse_complement = ''.join([complement[base] for base in reversed_sequence])
#   return reverse_complement

# fasta_file = 'lambda_virus.fa'

# pattern = 'AGTCGA'

# leftmost_occurrence = find_leftmost_occurrence_in_fasta(fasta_file, pattern)

# if leftmost_occurrence == -1:
#   pattern = reverse_complement(pattern)
#   leftmost_occurrence = find_leftmost_occurrence_in_fasta(fasta_file, pattern)

# print(leftmost_occurrence)




# To implement the naive_2mm() function, we can modify the original naive() function to allow up to two mismatches per occurrence. The following Python code shows the implementation of the naive_2mm() function:


import re

def naive_2mm(pattern, text):
  """Returns a list of indices where the pattern occurs in the text with up to 2 mismatches."""

  matches = []
  for i in range(len(text) - len(pattern) + 1):
    mismatches = 0
    for j in range(len(pattern)):
      if text[i + j] != pattern[j]:
        mismatches += 1
    if mismatches <= 2:
      matches.append(i)
  return matches

fasta_file = 'lambda_virus.fa'

pattern = 'AGGAGGTT'

matches = []
with open(fasta_file, 'r') as f:
  for line in f:
    if line.startswith('>'):
      continue
    else:
      matches.extend(naive_2mm(pattern, line))

print(len(matches))



import re

def count_occurrences(fasta_file, pattern):
  """Returns the total number of times the pattern occurs in the FASTA file."""

  total_occurrences = 0
  with open(fasta_file, 'r') as f:
    for line in f:
      if line.startswith('>'):
        continue
      else:
        occurrences = re.findall(pattern, line)
        total_occurrences += len(occurrences)
  return total_occurrences


def reverse_complement(sequence):
  """Returns the reverse complement of the sequence."""

  complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
  reversed_sequence = sequence[::-1]
  reverse_complement = ''.join([complement[base] for base in reversed_sequence])
  return reverse_complement
print(reverse_complement('AGTCGA'))
fasta_file = 'lambda_virus.fa'

pattern = 'AGGT'
pattern2 = 'ACCT'
# Count the number of times AGGT occurs in the FASTA file.
count_aggt = count_occurrences(fasta_file, pattern)
# Count the number of times ACCT occurs in the FASTA file.
count_acct = count_occurrences(fasta_file, pattern2)
# Get the total number of occurrences of AGGT and ACCT.
total_occurrences = count_aggt + count_acct

print(total_occurrences)





####################################

def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome
genome = readGenome('lambda_virus.fa')
print(genome[:1])

def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t

reverse=reverseComplement(genome)


def naive(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences
s=(len(naive('TTAA',genome)))
a=(len(naive('TTAA',reverse)))
print(s)
print(a)
print(a+s)




# occurrences = naive('AGTCGA', genome) 
# occurrences1 = naive('TCGACT', reverse)
# # Find the offset of the leftmost occurrence
# offset1 = min(occurrences)
# offset=min(occurrences1)

# # Print the offset
# print(offset)
# print(offset1)
def leftmost(p, t):
  occurences=naive(p,t)
  if occurences :
    return min(occurences)
  else:
    return None



# Read the Lambda virus genome
genome = readGenome('lambda_virus.fa')

# Find the leftmost occurrence of ACTAAGT or its reverse complement
def print_leftmost(p, t):
  
  leftmost_occurrence = leftmost(p, t)
  if leftmost_occurrence is not None:
    
    print('Leftmost occurrence of {}: {}'.format(p, leftmost_occurrence))
  else:
    print('{} not found in the Lambda virus genome'.format(p))

# Print the leftmost occurrence of ACTAAGT or its reverse complement
print_leftmost('TCGACT', genome)




def naive_2mm(p, t):
  """Finds all occurrences of a pattern in a text with up to 2 mismatches."""
  occurrences = []
  for i in range(len(t) - len(p) + 1):
    mismatches = 0
    for j in range(len(p)):
      if t[i+j] != p[j]:
        mismatches += 1
    if mismatches <= 2:
      occurrences.append(i)
  return occurrences

# Read the Lambda virus genome
genome = readGenome('lambda_virus.fa')

# # Find all occurrences of TTCAAGCC in the Lambda virus genome with up to 2 mismatches
occurrences = naive_2mm('TTCAAGCC', genome)

# # Print the number of occurrences
print(len(occurrences))

occurrences = naive_2mm('AGGAGGTT', genome) + naive_2mm('AACCTCCT', reverse)

# Find the offset of the leftmost occurrence
offset = min(occurrences)

# Print the offset
print(offset)
