#!/usr/bin/env python3
"""Independent Python exhaustive checks, synthetic only, not biological evidence."""
import csv, io, gzip, pathlib, random, subprocess, tempfile, unittest
class OracleTest(unittest.TestCase):
 def test_exhaustive_counts(self):
  rng=random.Random(20260925)
  seqs=[''.join(rng.choice('ACGT') for _ in range(19)) for _ in range(36)]
  seqs += [seqs[0][:-1]+next(b for b in 'ACGT' if b!=seqs[0][-1]),seqs[1][1:]+'A']
  with tempfile.TemporaryDirectory() as d:
   d=pathlib.Path(d);lib=d/'library.tsv';lib.write_text('id\tsequence\tgene\n'+''.join(f'g{i}\t{s}\tG{i//3}\n' for i,s in enumerate(seqs)))
   reads=[]
   for s in seqs:
    reads+=['T'*23+s+'AACG'*3,'T'*22+s+'AACG'*3]
    for k in (1,2):
     t=list(s)
     for j in rng.sample(range(19),k):t[j]=next(b for b in 'ACGT' if b!=t[j])
     reads.append('T'*23+''.join(t)+'AACG')
    reads.append('T'*23+s[:6]+'N'+s[7:]+'AACG')
   reads+=['N'*42,'ACGT','T'*42,'T'*23+seqs[0]+'AACGT']
   with gzip.open(d/'reads.fq.gz','wt') as f:
    for i,s in enumerate(reads):f.write(f'@{i}\n{s}\n+\n'+'I'*len(s)+'\n')
   binary=pathlib.Path(__file__).resolve().parent/'oracle'
   subprocess.run([str(binary),str(lib),str(d/'reads.fq.gz'),str(d/'out'),'100'],check=True)
   got=list(csv.DictReader(io.StringIO((d/'out.counts.tsv').read_text()),delimiter='\t'))
   offsets={mode:[int(r['offset']) for r in csv.DictReader(io.StringIO((d/'out.offsets.tsv').read_text()),delimiter='\t') if r['mode']==mode and r['retained']=='1'] for mode in ['exact','best']}
   expected={k:[0]*len(seqs) for k in list(got[0])[2:]}
   def hit(s,o,allow=True):
    w=s[o:o+19]
    if len(w)!=19 or (not allow and any(b not in 'ACGT' for b in w)):return [],[]
    ds=[sum(a!=b for a,b in zip(w,t)) for t in seqs]
    return [i for i,x in enumerate(ds) if x==0],[i for i,x in enumerate(ds) if x<=1]
   for s in reads:
    e,r=hit(s,23);b=e if e else r
    for k,ids in [('fixed_exact',e),('fixed_radius1',r),('fixed_best1',b)]:
     if len(ids)==1:expected[k][ids[0]]+=1
    for mode in ['exact','best']:
     all_ids=[];joint_e=set();joint_r=set()
     for o in offsets[mode]:
      e,r=hit(s,o,False);ids=e if mode=='exact' else e or r
      if len(ids)==1:all_ids.extend(ids)
      joint_e.update(e);joint_r.update(r)
     suffix='exact' if mode=='exact' else 'best1'
     for i in all_ids:expected['multi_'+suffix][i]+=1
     for i in set(all_ids):expected['distinct_'+suffix][i]+=1
     joint=joint_e or joint_r
     if mode=='best' and len(joint)==1:expected['joint_best1'][next(iter(joint))]+=1
   for k,v in expected.items():self.assertEqual([int(r[k]) for r in got],v,k)
   states=list(csv.DictReader(io.StringIO((d/'out.qc.tsv').read_text()),delimiter='\t'))
   for mode in ['fixed_exact','fixed_radius1','fixed_best1']:self.assertEqual(sum(int(r['value']) for r in states if r['mode']==mode),len(reads))
   self.assertEqual(sum(expected['fixed_exact']),sum(len(hit(s,23)[0])==1 for s in reads))
   print('Validated eight policy matrices and three state denominators on',len(reads),'synthetic reads.')
if __name__=='__main__':unittest.main()
