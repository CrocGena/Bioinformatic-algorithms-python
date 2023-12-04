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
print(edge)

def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's suffx in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match

import itertools
def scs(ss):
    """ Returns shortest common superstring of given
        strings, which must be the same length """
    shortest_sup = None
    for ssperm in itertools.permutations(ss):
        sup = ssperm[0]  # superstring starts as first string
        for i in range(len(ss)-1):
            # overlap adjacent strings A and B in the permutation
            olen = overlap(ssperm[i], ssperm[i+1], min_length=1)
            # add non-overlapping portion of B to superstring
            sup += ssperm[i+1][olen:]
        if shortest_sup is None or len(sup) < len(shortest_sup):
            shortest_sup = sup  # found shorter superstring
    return shortest_sup  # return shortest

def scs_list(ss):
    """ Returns the alphebeticaly sorted list of shortest common superstrings of given strings, which must be the same length """
    import itertools
    shortest_sup = []
    shortest_length = 0
    for ssperm in itertools.permutations(ss):
        sup = ssperm[0]  # superstring starts as first string
        for i in range(len(ss)-1):
            # overlap adjacent strings A and B in the permutation
            olen = overlap(ssperm[i], ssperm[i+1], min_length=1)
            # add non-overlapping portion of B to superstring
            sup += ssperm[i+1][olen:]
        if shortest_length == 0:
            shortest_sup.append(sup)
            shortest_length = len(sup)
        elif len(sup) < shortest_length:
            shortest_length = len(sup)
            shortest_sup = [sup]
        elif len(sup) == shortest_length:
            shortest_sup.append(sup)
        else:
            #simply ignore and move on
            None
    return sorted(shortest_sup)

strin=['CCT','CTT','TGC','TGG','GAT','ATT']
shortle=scs(strin)
print(shortle)
print(len(shortle))


shortest=scs_list(strin)
print(len(shortest))