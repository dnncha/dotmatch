import csv,gzip,json,os,pathlib,random,subprocess,sys,tempfile,unittest
P=pathlib.Path(__file__).resolve().parent
sys.path.insert(0,str(P.parent))
import joint_resolution as ar3
BIN=os.environ.get('AR004_COUNTER',str(P/'count'))
class Tests(unittest.TestCase):
 def test_all_position_oracle(self):
  rng=random.Random(20260906);comparisons=0
  with tempfile.TemporaryDirectory() as tmp:
   p=pathlib.Path(tmp)/'library.tsv'
   for rep in range(16):
    ss=[''.join(rng.choices('ACGTN',k=4)) for _ in range(10)];ss[1]=ss[0]
    rows=[dict(id=str(i),sequence=s,gene=str(i//2)) for i,s in enumerate(ss)]
    with p.open('w') as f:
     w=csv.DictWriter(f,fieldnames=['id','sequence','gene'],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    seqs=[''.join(rng.choices('ACGTN',k=6)) for _ in range(500)]+['ACGT','ACGT??','AAAAAA']
    dec=ar3.Decoder(rows,[0,1,2],cache_size=0)
    r=subprocess.run([BIN,str(p),'0,1,2','--probe'],input='\n'.join(seqs)+'\n',text=True,capture_output=True,check=True)
    lines=r.stdout.splitlines();self.assertEqual(len(lines),len(seqs))
    for seq,line in zip(seqs,lines):
     pieces=line.split('\t');valid=int(pieces[0]);hits=tuple(tuple(map(int,x.split(','))) for x in pieces[1].split(';') if x)
     expected=ar3.exhaustive(dec,seq);self.assertEqual(valid,expected is not None)
     self.assertEqual(hits,() if expected is None else expected)
     calls=[tuple(int(x) for x in v.split(',') if x) for v in pieces[2:]]
     _,oracle=dec.decode(seq)
     self.assertEqual(calls[2:],[() if oracle[i] is None else oracle[i][0] for i in (0,2,1)])
     if expected is not None:
      event=[[],[]]
      for off in (0,1,2):
       h=[(i,d) for o,i,d in expected if o==off];e=[i for i,d in h if d==0];bd=min((v for _,v in h),default=2);b=[i for i,d in h if d==bd]
       if len(e)==1:event[0]+=e
       if len(b)==1:event[1]+=b
      self.assertEqual(calls[:2],[tuple(x) for x in event])
     comparisons+=1
  print('Exhaustive constructed all-target/all-position comparisons:',comparisons)
 def test_stream_cache_thinning_and_failure(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=pathlib.Path(tmp);lib=p/'l.tsv';lib.write_text('id\tsequence\tgene\na\tACGT\tg\nb\tCGTA\tg\n');reads=p/'r.fastq.gz'
   with gzip.open(reads,'wt') as f:
    for i in range(100):f.write(f'@{i}\nACGTA\n+\nIIIII\n')
   for cache in (0,1,16):subprocess.run([BIN,str(lib),'0,1',str(reads),str(p/str(cache)),'.5','42',str(cache)],check=True,capture_output=True)
   for name in ('guide-counts.tsv','gene-counts.tsv','qc.tsv','metrics.json','shared-read-pairs.tsv'):
    self.assertEqual((p/'0'/name).read_bytes(),(p/'16'/name).read_bytes());self.assertEqual((p/'0'/name).read_bytes(),(p/'1'/name).read_bytes())
   self.assertEqual(json.loads((p/'0/metrics.json').read_text())['policies']['event_best']['extra_events'],100)
   with gzip.open(p/'broken.gz','wt') as f:f.write('@x\nACGT\n+\n')
   self.assertNotEqual(subprocess.run([BIN,str(lib),'0,1',str(p/'broken.gz'),str(p/'bad'),'1','42','1'],capture_output=True).returncode,0)
   self.assertNotEqual(subprocess.run([BIN,str(lib),'0,1',str(reads),str(p/'0'),'1','42','1'],capture_output=True).returncode,0)
if __name__=='__main__':unittest.main(verbosity=2)
