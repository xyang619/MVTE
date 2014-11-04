'''
Name: macssim2eigen.py
Date: 2014-9-9
Last modified: 2014-09-15
Version: 1.02
Author: Young
Description:
    Convert the output of MACS simulation into eigenstrat format
Input file: the standard output (stdout) of MACS simulation
Output files: prefix.ind prefix.snp prefix.geno
Arguments:
    -h --help               print help
    -f --file filename      name of MACS standard output file [string]
    -n --ninds n1,n2 ...    #samples in each subpopulation [integer]
    -l --length L           length of simulated sequence [integer]
    -p --prefix prefix      prefix of output files [string]
'''

import sys, random, getopt

def gen_allele():
    return random.sample('AGCT',2)

def gen_ind(ninds, prefix):
    indf = '{}.ind'.format(prefix)
    with open(indf, 'w') as out:
        i = 1
        for nind in ninds:
            for k in range(nind):
                out.write('SAM{}_{}\tU\tPOP{}\n'.format(i,k+1,i))
            i+=1
    print('Write indfile into {}.ind'.format(prefix))

#SITE:	0	   0.0279543585	     0.17045867	1001000100
def gen(simfile, ninds=None ,L=1e7, prefix='sim'):
    gen_ind(ninds, prefix)
    snpf = '{}.snp'.format(prefix)
    genf = '{}.geno'.format(prefix)
    with open(simfile) as f, open(snpf, 'w') as sf, open(genf, 'w') as gf:
        chrom=1
        for line in f:
            if line.startswith('SITE'):
                lst=line.split()
                rsid='rs'+lst[1]
				pp=int(float(lst[2])*L)
                gp=float(lst[2])*L/100000000.0	#convert to Morgan
                gp='{:.8f}'.format(gp)
                a1,a2=gen_allele()
                sf.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(rsid,chrom,gp,pp,a1,a2))
                gn=lst[4] 
                gl=''
                for i in range(len(gn)/2):
                    gl+=repr(int(gn[2*i])+int(gn[2*i+1]))
                gf.write(gl+'\n')
            elif line.startswith('END_SELECTED_SITES'):
                chrom+=1
    print('Write snpfile into {}.snp'.format(prefix))
    print('write genofile into {}.geno'.format(prefix))

def usage():
    print('''Description:
    Convert the output of MACS simulation into eigenstrat format
Input file: the standard output (stdout) of MACS simulation
Output files: prefix.ind prefix.snp prefix.geno
Arguments:
    -h --help               print help
    -f --file filename      name of MACS standard output file [string]
    -n --ninds n1,n2 ...    #samples in each subpopulation [integer]
    -l --length L           length of simulated sequence [integer]
    -p --prefix prefix      prefix of output files [string]
''')

def main():
    f=''
    ns=[]
    l=1e7
    p='sim'
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hf:n:l:p:', 
        ['help','file','npops','length','prefix'])
    except getopt.GetoptError as err:
        print(err)
        usage()
        sys.exit(2)
    for o, a in opts:
        if o in ('-h', '--help'):
            usage()
            sys.exit(1)
        elif o in ('-f', '--file'):
            f=a
        elif o in ('-n', '--ninds'):
            for k in a.split(','):
                ns.append(int(k))
        elif o in ('-l', '--length'):
            l=int(a)
        elif o in ('-p', '--prefix'):
            p=a
    #print(f,n,ns,l,p)
    assert (len(f) > 0), 'Input file is empty'
    gen(f, ns, l, p)    

if __name__=='__main__':
    main()
                              
