"""peca summary table and plots"""
from collections import Counter
from pathlib import Path
from sys import stdout

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


NPROT = sum(1 for line in open("x.txt"))-1
NSAMPLES = sum(1 for line in open("s_RR"))


NREPS, NTIME = (int(x) for x in open("attr.txt").readline().rstrip().split('\t')[:2])
#print(NREPS, NTIME)


RR = np.zeros((NPROT, NTIME-1))
for line in (l.rstrip().split('\t') for l in open("s_RR")):
    RR += np.array(line, dtype=float).reshape(NPROT, NTIME-1)
RR /= NSAMPLES


CPR = np.zeros((NPROT, NTIME-2))
for line in (l.rstrip().split('\t') for l in open("s_CPR")):
    CPR += np.array(line, dtype=int).reshape(NPROT, NTIME-2)
CPR /= NSAMPLES


DD = np.zeros((NPROT, NTIME-1))
for line in (l.rstrip().split('\t') for l in open("s_DD")):
    DD += np.array(line, dtype=float).reshape(NPROT, NTIME-1)
DD /= NSAMPLES


CPD = np.zeros((NPROT, NTIME-2))
for line in (l.rstrip().split('\t') for l in open("s_CPD")):
    CPD += np.array(line, dtype=int).reshape(NPROT, NTIME-2)
CPD /= NSAMPLES

### append the fdr columns
#for syn
CNT = Counter(CPR.flatten())

SORTKEY = sorted(CNT, reverse=True)

FDRMAP = dict()
DENOM = CNT[SORTKEY[0]]
NUMER = (1-SORTKEY[0])*DENOM
for key in SORTKEY:
    FDRMAP[key] = NUMER/DENOM
    NUMER += (1-key)*CNT[key]
    DENOM += CNT[key]

#for deg
CNT = Counter(CPD.flatten())

SORTKEY = sorted(CNT, reverse=True)

FDRMAP_D = dict()
DENOM = CNT[SORTKEY[0]]
NUMER = (1-SORTKEY[0])*DENOM
for key in SORTKEY:
    FDRMAP_D[key] = NUMER/DENOM
    NUMER += (1-key)*CNT[key]
    DENOM += CNT[key]

MPATH = Path("M.txt").is_file()

with open('normX.txt') as xX, \
        open('normH.txt') as hH, \
        (open('normM.txt') if MPATH else open('normH.txt')) as mM, \
        open('data_R_CPS.txt', 'w') as cp_:
    #cp_.write(xX.readline().rstrip().split('\t', 1)[1] \
    #        +('\t'+mM.readline().rstrip().split('\t', 1)[1] if MPATH else '') \
    #        +'\t'+hH.readline().rstrip().split('\t', 1)[1] \
    cp_.write('\t'.join(z+"_X"+str(int(n/NTIME))+"t"+str(n%NTIME) \
        for n, z in enumerate(xX.readline().rstrip().split('\t')[1:])) \
        +('\t'+'\t'.join(z+"_M"+str(int(n/NTIME))+"t"+str(n%NTIME) \
        for n, z in enumerate(mM.readline().rstrip().split('\t')[1:])) if MPATH else '') \
        +'\t'+'\t'.join(z+"_Y"+str(int(n/NTIME))+"t"+str(n%NTIME) \
        for n, z in enumerate(hH.readline().rstrip().split('\t')[1:])) \
        +'\t'+'\t'.join(['R'+str(i) for i in range(NTIME-1)]) \
        +'\t'+'\t'.join(['D'+str(i) for i in range(NTIME-1)]) \
        +'\t'+'\t'.join(['signedCPS'+str(i) for i in range(1, NTIME-1)]) \
        +'\t'+'\t'.join(['signedCPD'+str(i) for i in range(1, NTIME-1)]) \
        +'\t'+'\t'.join(['FDR_S'+str(i) for i in range(1, NTIME-1)]) \
        +'\t'+'\t'.join(['FDR_D'+str(i) for i in range(1, NTIME-1)]) \
        +'\n')
    for i in range(NPROT):
        cp_.write(xX.readline().rstrip() \
            +('\t'+mM.readline().rstrip().split('\t', 1)[1] if MPATH else '') \
            +'\t'+hH.readline().rstrip().split('\t', 1)[1] \
            +'\t'+'\t'.join(str(x) for x in RR[i,]) \
            +'\t'+'\t'.join(str(x) for x in DD[i,]) \
            +'\t'+'\t'.join(str(x*(1 if RR[i, n] > RR[i, n-1] else -1)) \
            for n, x in enumerate(CPR[i,], 1)) \
            +'\t'+'\t'.join(str(x*(1 if RR[i, n] > RR[i, n-1] else -1)) \
            for n, x in enumerate(CPD[i,], 1)) \
            +'\t'+'\t'.join(str(FDRMAP[x]) for x in CPR[i,]) \
            +'\t'+'\t'.join(str(FDRMAP_D[x]) for x in CPD[i,]) \
            +'\n')


###loglikelihood traceplot########################
#plt.figure().set_size_inches(9, 9)
plt.title('loglikelihood traceplot')
plt.plot(range(1, NSAMPLES+1), [float(x) for x in open('s_loglike').readlines()])
plt.savefig('trace_loglike.pdf')
##################################################


EFSH = dict()
XAX = []
with open('EfsH.txt') as efsh:
    XAX = efsh.readline().rstrip('\n').split('\t')[3:]
    for line in (l.rstrip('\n').split('\t') for l in efsh):
        EFSH['\t'.join(line[:3])] = line[3:]

EFSM = dict()
if MPATH:
    with open('EfsM.txt') as efsm:
        efsm.readline()
        for line in (l.rstrip('\n').split('\t') for l in efsm):
            EFSM['\t'.join(line[:3])] = line[3:]

EFSX = dict()
EFSXPATH = Path("EfsX.txt").is_file()
if EFSXPATH:
    with open('EfsX.txt') as efsx:
        efsx.readline()
        for line in (l.rstrip('\n').split('\t') for l in efsx):
            EFSX['\t'.join(line[:3])] = line[3:]

ETAH = np.array(open('mean_etaH').readline().rstrip().split('\t'), \
        dtype=float).reshape(NPROT, NTIME*NREPS)

if MPATH:
    ETAM = np.array(open('mean_etaM').readline().rstrip().split('\t'), \
            dtype=float).reshape(NPROT, NTIME*NREPS)


with PdfPages('mRNAprot.pdf') as pdf, \
        open('X.txt') as xx, \
        open('H.txt') as hh, \
        (open('M.txt') if MPATH else open('H.txt')) as mm, \
        open('normX.txt') as xX, \
        open('normH.txt') as hH, \
        (open('normM.txt') if MPATH else open('normH.txt')) as mM:
    PR = 2
    LOGHY = 'log(protein) '
    if MPATH:
        PR = 3
        mm.readline()
        mM.readline()
        LOGHY = 'log(H) '
    xx.readline()
    hh.readline()
    xX.readline()
    hH.readline()
    for p, linex in enumerate(xx):
        x = [(i if i != 'NA' else np.nan) for i in linex.rstrip().split('\t')[1:]]
        x = np.array(x, dtype=float)
        X = np.array(xX.readline().rstrip().split('\t')[1:], dtype=float)
        h = [(i if i != 'NA' else np.nan) for i in hh.readline().rstrip().split('\t')[1:]]
        h = np.array(h, dtype=float)
        H = np.array(hH.readline().rstrip().split('\t')[1:], dtype=float)
        if MPATH:
            m = [(i if i != 'NA' else np.nan) for i in mm.readline().rstrip().split('\t')[1:]]
            m = np.array(m, dtype=float)
            M = np.array(mM.readline().rstrip().split('\t')[1:], dtype=float)
        if p > 10:
            break
        plt.figure(figsize=((NREPS+2)*4, PR*4))
        plt.suptitle(linex.rstrip().split('\t', 1)[0])
        print('\x08'*99, p, '/', NPROT, end=' ')
        stdout.flush()

        for j in range(NREPS):
            plt.subplot(PR, NREPS+2, j+1).set_title('log(mRNA) '+str(j+1))
            plt.ylim(np.nanmin(x), np.nanmax(x))
            plt.scatter(range(NTIME), x[NTIME*j:NTIME*(j+1)], c='k')
            plt.scatter(range(NTIME), X[NTIME*j:NTIME*(j+1)], \
                    facecolors='none', \
                    edgecolors=np.where(np.isnan(x[NTIME*j:NTIME*(j+1)]), 'red', 'black'))
            if EFSXPATH:
                plt.plot(XAX, EFSX[str(p)+'\t0\t'+str(j)], 'k-')

        plt.subplot(PR, NREPS+2, NREPS+1).set_title('Protein synthesis')
        plt.xlim(0, NTIME-1)
        plt.step(np.arange(0.5, NTIME-0.5), RR[p,], 'k-', where='mid')
        plt.ylim(0)

        plt.subplot(PR, NREPS+2, NREPS+2).set_title('Protein degradation')
        plt.xlim(0, NTIME-1)
        plt.step(np.arange(0.5, NTIME-0.5), DD[p,], 'k-', where='mid')
        plt.ylim(0)

        for j in range(NREPS):
            plt.subplot(PR, NREPS+2, NREPS+3+j).set_title(LOGHY+str(j+1))
            plt.ylim(np.nanmin(h), np.nanmax(h))
            plt.scatter(range(NTIME), h[NTIME*j:NTIME*(j+1)], c='k')
            plt.scatter(range(NTIME), H[NTIME*j:NTIME*(j+1)], \
                    facecolors='none', \
                    edgecolors=np.where(np.isnan(h[NTIME*j:NTIME*(j+1)]), 'red', 'black'))
            plt.plot(range(NTIME), ETAH[p, NTIME*j:NTIME*(j+1)], 'k--')
            plt.plot(XAX, EFSH[str(p)+'\t0\t'+str(j)], 'k-')

        plt.subplot(PR, NREPS+2, 2*NREPS+3).set_title('Change point probability')
        plt.xlim(0, NTIME-1)
        plt.ylim(0, 1)
        objects = range(1, NTIME-1)
        y_pos = np.arange(1, len(objects)+1)
        plt.bar(y_pos, CPR[p,], align='center', alpha=0.5)
        plt.xticks(y_pos, objects)

        plt.subplot(PR, NREPS+2, 2*NREPS+4).set_title('Change point probability')
        plt.xlim(0, NTIME-1)
        plt.ylim(0, 1)
        plt.bar(y_pos, CPD[p,], align='center', alpha=0.5)
        plt.xticks(y_pos, objects)

        if MPATH:
            for j in range(NREPS):
                plt.subplot(PR, NREPS+2, 2*NREPS+5+j).set_title('log(M) '+str(j+1))
                plt.ylim(np.nanmin(m), np.nanmax(m))
                plt.scatter(range(NTIME), m[NTIME*j:NTIME*(j+1)], c='k')
                plt.scatter(range(NTIME), M[NTIME*j:NTIME*(j+1)], \
                        facecolors='none', \
                        edgecolors=np.where(np.isnan(m[NTIME*j:NTIME*(j+1)]), 'red', 'black'))
                plt.plot(range(NTIME), ETAM[p, NTIME*j:NTIME*(j+1)], 'k--')
                plt.plot(XAX, EFSM[str(p)+'\t0\t'+str(j)], 'k-')

        pdf.savefig()  # saves the current figure into a pdf page
        plt.close()

    D = pdf.infodict()
    D['Title'] = 'gene-specific plots'
    D['Author'] = 'guo shou TEO'
