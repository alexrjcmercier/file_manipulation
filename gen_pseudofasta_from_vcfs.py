#coding:utf-8
#!/bin/env python2.7
#All inputs and outputs are hardcoded into the script
#It uses the individual/vcf file names and the chromosome number to iterate on vcf files and find the line corresponding to the different chromosomes.
#It also need a reference genome file

vcffiles = ['G1_2', 'G2', 'G3', 'G4', 'G5', 'G6', 'G7', 'G8', 'G9', 'G10', 'G11', 'G12', 'G13', 'G14', 'G15', 'G16', 'H1', 'R1', 'R2', 'T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'T8', 'T9', 'T10', 'T11', 'T12', 'T13', 'Bpse', 'Bfab']
chrlist = ['BCIN01', 'BCIN02', 'BCIN03', 'BCIN04', 'BCIN05', 'BCIN06', 'BCIN07', 'BCIN08', 'BCIN09', 'BCIN10', 'BCIN11', 'BCIN12', 'BCIN13', 'BCIN14', 'BCIN15', 'BCIN16', 'BCIN17', 'BCIN18']

referencefile = '/work/daphne/24genomes2017/BcinB0510_finalassembly_January2015.fasta'
inputdir = '/work/daphne/24genomes2017/3-snpcalling/selected-edit/'
outputdir = '/work/daphne/24genomes2017/5-sweeps/genovar_stats/seqID/'

fastadict = dict()
with open (referencefile, 'r') as reffile :
    fastaline = reffile.readline()
    while fastaline != '' :
        if '>' in fastaline :
            chrID = fastaline[1:7]
            fastadict[chrID] = ''
        else:
            fastadict[chrID] = ''.join([str(x) for x in [fastadict[chrID], (fastaline.replace('\n', ''))]])
        fastaline = reffile.readline()

sizedict = dict()
for k, v in sorted(fastadict.items()) :
    sizedict[k] = len(v)
    with open(''.join([str(x) for x in [outputdir, 'seqfile_', k, '.fasta']]), 'w') as outputfile :
        refbuffer = ''.join([str(x) for x in ['>', 'B05.10', '\n', v, '\n']])
        outputfile.write(refbuffer)

for chrid in chrlist :
    with open(''.join([str(x) for x in [outputdir, 'seqfile_', chrid, '.fasta']]), 'a') as outputfile :
        for vcf in vcffiles :
            with open(''.join([str(x) for x in [inputdir, 'selected-', vcf, '-edit.vcf']]), 'r') as vcffile :
                fastabuffer = str()
                posindex = 1
                linebuffer = vcffile.readline().split('\t')
                while '#' in linebuffer[0] : #skip header lines
                    linebuffer = vcffile.readline().split('\t')
                while str(linebuffer[0]) != str(chrid) : #skip wrong chromosome lines
                    linebuffer = vcffile.readline().split('\t')
                for nucpos in range(1, int(sizedict[chrid]+1), 1) : #1;n+1 because python iterates on 0;n-1, whereas nucleotide position needs to be 1;n
                    if chrid in linebuffer :
                        if nucpos != int(linebuffer[1]) :
                            fastabuffer += (str(fastadict[chrid][nucpos-1]))
                        elif nucpos == int(linebuffer[1]) :
                            fastabuffer += (str(linebuffer[4]))
                            linebuffer = vcffile.readline().split('\t')
                    else :
                        fastabuffer += (str(fastadict[chrid][nucpos-1]))
                vcfbuffer = ''.join([str(x) for x in ['>', vcf, '\n', fastabuffer, '\n']])
                outputfile.write(vcfbuffer)
