import gzip, sys, numpy as np
from scipy.stats import spearmanr
def load(p, col):
    f=gzip.open(p,'rt'); h=f.readline().rstrip('\n').split('\t'); i=h.index(col); iq=h.index('n_tx_quant') if 'n_tx_quant' in h else None
    return {x[0]:(float(x[i]), int(x[iq]) if iq else 1) for x in (l.rstrip('\n').split('\t') for l in f)}
kr=load(sys.argv[1],'est_counts'); rs=load(sys.argv[2],'est_reads')
g=[k for k in kr if kr[k][0]>=50 and k in rs and rs[k][1]>0]
a=np.array([kr[k][0] for k in g]); b=np.array([rs[k][0] for k in g])
print('genes kallisto>=50:',sum(v[0]>=50 for v in kr.values()),' with quantifiable rnaskim tx:',len(g),' not in rnaskim index:',sum(1 for k in kr if kr[k][0]>=50 and k not in rs), ' unquantifiable in index:',sum(1 for k in kr if kr[k][0]>=50 and k in rs and rs[k][1]==0))
print('Spearman',spearmanr(a,b)[0],' total kallisto',sum(v[0] for v in kr.values()),' total rnaskim',sum(v[0] for v in rs.values()))
