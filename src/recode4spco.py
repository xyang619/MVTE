'''Recode the genotype in eigenstrat format into StepPCO format'''

import sys

def code(allele):
    if allele=='0':
        return '1'
    elif allele=='1':
        return '0'
    elif allele=='2':
        return '-1'
    else:
        return 'NA'

def convert(infile,outfile):
    with open(infile) as f, open(outfile,'w') as out:
        for line in f:
            outStr=''
            for gp in line[:-1]:
                outStr+=code(gp)+'\t'
            outStr+='\n'
            out.write(outStr)
    print('Finished converting')
if __name__=='__main__':
    if len(sys.argv)>2:
        convert(sys.argv[1],sys.argv[2])
    else:
        print('Usage: python {} infile outfile'.format(sys.argv[0]))
        
    
