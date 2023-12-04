def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome


sequences =readGenome('chr1.GRCh38.excerpt.fasta')
print(sequences[:100])





def eddist(a,b):
    if  len(a)==0:
        return len(b)
    if len(b)==0:
        return len(a)
    delt=1 if a[-1]!= b[-1] else 0
    return min(eddist(a[:-1],b[:-1])+delt , eddist(a,b[:-1])+1,eddist(a[:-1],b)+1)
# import datetime as d
# st=d.datetime.now()
# eddist("Shakespeare","sake spear")
# stt=d.datetime.now()-st
# total_seconds =stt.total_seconds()
# print(total_seconds)
    
    
def edistDistance(x,y):
    D=[] # matrix create where fill in edit distances substrings where already computed
    for i in range(len(x)+1):
        D.append([0]*(len(y)+1))
    for i in range(len(x)+1):
        D[i][0]=i
    
    for i in range(len(y)+1):
        D[0][i]=i
    for i in range(1,len(x)+1):
        for j in range(1,len(y)+1):
            disthor=D[i][j-1]+1
            distver=D[i-1][j]+1
            if x[i-1]== y[j-1]:
                distdiag=D[i-1][j-1]
            else:
                distdiag=D[i-1][j-1]+1
            D[i][j]=min(disthor,distver,distdiag)
    return D[-1][-1]


import timeit


x = 'GCTGATCGATCGTACG'
y = sequences
edit_distance = edistDistance(x, y)





execution_time = timeit.timeit(lambda: edistDistance(x, y), number=1)

# Print the edit distance and the time taken.
print(f"Edit distance: {edit_distance}")
print(f"Time taken: {execution_time} seconds")

# global aligment function with penaly matrix
alphabet=['A','C','G','T']
Score=[[0,4,2,4,8],
       [4,0,4,2,8],
       [2,4,0,4,8],
       [4,2,4,0,8],
       [8,8,8,8,8]]

def GlobalAlignment(x,y):
    D=[] # matrix create where fill in edit distances substrings where already computed
    for i in range(len(x)+1):
        D.append([0]*(len(y)+1))
    for i in range(1,len(x)+1):
        D[i][0]=D[i-1][0]+Score[alphabet.index(x[i-1])][-1]
    
    for i in range(1,len(y)+1):
        D[0][i]=D[0][i-1]+Score[-1][alphabet.index(y[i-1])]
    for i in range(1,len(x)+1):
        for j in range(1,len(y)+1):
            disthor=D[i][j-1]+Score[-1][alphabet.index(y[j-1])]
            distver=D[i-1][j]+Score[alphabet.index(x[i-1])][-1]
            if x[i-1]== y[j-1]:
                distdiag=D[i-1][j-1]
            else:
                distdiag=D[i-1][j-1]+Score[alphabet.index(x[i-1])][alphabet.index(y[j-1])]
            D[i][j]=min(disthor,distver,distdiag)
    return D[-1][-1]


x = 'AAGCTTAA'
y = 'AAGCT'
GlobalAlignments = GlobalAlignment(x, y)
execution_time = timeit.timeit(lambda: GlobalAlignment(x, y), number=1)              

print(f"\nGlobalAlignments: {GlobalAlignments}")
print(f"Time taken: {execution_time} seconds")          
            
        
def overlap(a, b, min_len=3):
    start = 0
    while True:
        start = a.find(b[:min_len], start)
        if start == -1:
            return 0
        if b.startswith(a[start:]):
            return len(a) - start
        start += 1
overla=overlap('sadafdamdaCGT','CGTACCGGT')
print(overla)
       
       

from itertools import permutations

# print(list(permutations([1,2,3],1)))


def naice_overlap_map(reads,k):
    olaps={}
    for a,b in permutations(reads,2):
        overlap_length =overlap(a,b,min_len=k)
        if overlap_length>0:
            olaps[(a,b)]=overlap_length
    return olaps

reads=['ACGGGATGATC','ATCAAGTGGA','GGAAAGTACGGA']
print(naice_overlap_map(reads,3))


def de_bruijin_ize(st,k):
    edges=[]
    node=set()
    for i in range(len(st)-k+1):
        edges.append((st[i:i+k-1],st[i+1:i+k]))
        node.add(st[i:i+k-1])
        node.add(st[i+1:i+k])
    return node,edges

node, edge=de_bruijin_ize('ACGCGTCG',3)

print(node)

        
