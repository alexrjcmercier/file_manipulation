import os, sys
scriptname, inputfile, outputpath, windowsize, endofseq = sys.argv

inputfile_name = inputfile.split('/')[-1]
chrsizedict = {'BCIN01':4109373, 'BCIN02':3341473, 'BCIN03':3226611, 'BCIN04':2468882, 'BCIN05':2959378, 'BCIN06':2725906, 
'BCIN07':2652353, 'BCIN08':2617329, 'BCIN09':2547566, 'BCIN10':2419276, 'BCIN11':2359939, 'BCIN12':2352958, 
'BCIN13':2257609, 'BCIN14':2138025, 'BCIN15':2027721, 'BCIN16':1969743}

#endofseq :'truncated' or 'filled': 'truncated' will ommit the last part of the sequence which did not fit into a window.
#'filled' will add monomorphic non-folded sites at the end of the sequence to fill the last window and allow calculation.
#Default behavior is truncated.
outgroup = 'Bfab'

windowsize = int(windowsize)
truncatedict = dict()
for key, value in chrsizedict.items() :
    if endofseq == 'filled' :
        truncval = value - value%windowsize + windowsize
    else :
        truncval = value - value%windowsize
    truncatedict[key] = truncval
chrsizedict = truncatedict

posindex = 1
with open(inputfile, 'r') as infile :
    with open(''.join(str(x) for x in [outputpath, inputfile_name[:-4], '_', endofseq, '_', ''.join([str(windowsize/1000), 'kb']),'.sf']), 'w') as outfile :
        outbuffer = '\t'.join(str(x) for x in ['position', 'x', 'n', 'folded\n'])
        outfile.write(outbuffer)
        inline = infile.readline()

        while '#' in inline :
            indinline = inline
            inline = infile.readline().strip()

        intup = tuple(inline.split('\t'))
        indintup = tuple(indinline.split('\t')[9::])
        chrnum, pos, refnuc, nucintup = intup[0], int(intup[1]), intup[3], intup[9::]
        altnuc = intup[4].split(',')
        maxposindex = int(truncatedict[chrnum])
        nbindividuals = len(indintup)-1

        while inline != '' and pos < maxposindex :
            inddict = dict()
            intup = tuple(inline.split('\t'))
            indintup = tuple(indinline.split('\t')[9::])
            chrnum, pos, refnuc, nucintup = intup[0], int(intup[1]), intup[3], intup[9::]
            altnuc = intup[4].split(',')

            if pos < maxposindex :
                while posindex < pos and posindex < maxposindex :
                    outbuffer = '\t'.join(str(x) for x in [posindex, '0', nbindividuals, '0\n'])
                    outfile.write(outbuffer)
                    posindex += 1

                for i in range(0, len(indintup)) :
                    if nucintup[i] == './.' :
                        inddict[indintup[i].split('.')[0]] = refnuc
                    else :
                        for j in range(0, len(altnuc)) :
                            if nucintup[i][0] == str(j+1) :
                                inddict[indintup[i].split('.')[0]] = altnuc[j]

                fold = None
                allelelist = list()
                alleledict = dict()
                for ind, nuc in inddict.items() :
                    if ind == outgroup :
                        ancestral_state = nuc
                    else :
                        allelelist.append(nuc)
                nout = len(allelelist)

                for allele in set(allelelist) :
                    alleledict[allele] = allelelist.count(allele)
                
                xout = 0
                if ancestral_state in allelelist :
                    fold = 0
                    for allele, occurences in alleledict.items() :
                        if allele != ancestral_state : 
                            xout += occurences
                else:
                    fold = 1
                    for allele, occurences in alleledict.items() :
                        if allele != ancestral_state : 
                            xout = occurences
                outbuffer = ''.join(str(x) for x in [pos, '\t', xout, '\t', nout, '\t', fold, '\n'])
                outfile.write(outbuffer)
                inline = infile.readline().strip()
                posindex += 1

        while posindex <= maxposindex :
            outbuffer = '\t'.join(str(x) for x in [posindex, '0', nbindividuals, '0\n'])
            outfile.write(outbuffer)
            posindex += 1
        

