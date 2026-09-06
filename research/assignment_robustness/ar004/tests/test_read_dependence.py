import itertools, json, math, tempfile, unittest, sys
from pathlib import Path
import numpy as np
from scipy import stats
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import read_dependence as m

def fixture(root, cross=False, reuse=True):
    rows=[dict(id='a',gene='A',upstream_count_events=100),dict(id='b',gene='B' if cross else 'A',upstream_count_events=100)]
    m.tsv(root/'all-guide-counts.tsv',rows)
    multi=[dict(target_ids=json.dumps(['a','b']),genes=json.dumps(['A','B'] if cross else ['A']),read_count=100)] if reuse else []
    m.tsv(root/'multiple-guide-read-classes.tsv',multi,['target_ids','genes','read_count'])
    s=dict(reads=200,matched_reads=100 if reuse else 200,count_events=200,extra_events=100 if reuse else 0,repeated_same_target_events=0,multiple_distinct_guides=100 if reuse else 0)
    manifest=dict(completion='complete',statistics=s,files={p.name:m.sha(p) for p in root.glob('*.tsv')})
    m.dump(root/'completion.json',manifest)

class Tests(unittest.TestCase):
    def test_identical_count_table_distinct_variance(self):
        with tempfile.TemporaryDirectory() as a,tempfile.TemporaryDirectory() as b:
            a=Path(a);b=Path(b);fixture(a,reuse=True);fixture(b,reuse=False)
            self.assertEqual((a/'all-guide-counts.tsv').read_bytes(),(b/'all-guide-counts.tsv').read_bytes())
            x=m.load_classes(a);y=m.load_classes(b)
            self.assertEqual(list(x['C']),list(y['C']));self.assertEqual(x['Q'][0],400);self.assertEqual(y['Q'][0],200)
            self.assertEqual(x['keff'][0],1);self.assertEqual(y['keff'][0],2)
    def test_cross_gene_dependencies_not_within_gene_inflation(self):
        with tempfile.TemporaryDirectory() as a:
            a=Path(a);fixture(a,cross=True);x=m.load_classes(a)
            self.assertTrue(np.array_equal(x['C'],x['Q']))
            covariance=x['B'].T@m.sparse.diags(x['ns'])@x['B']
            self.assertEqual(covariance.toarray().tolist(),[[100,100],[100,100]])
    def test_convolution_exhaustive(self):
        hist={1:3,2:2,3:1};law=m.exact_law(hist);counter={}
        weights=[1]*3+[2]*2+[3]
        for bits in itertools.product((0,1),repeat=6):
            d=sum(w*(2*b-1) for w,b in zip(weights,bits));counter[d]=counter.get(d,0)+1
        for d,p in zip(law['d'],law['pmf']):self.assertAlmostEqual(p,counter.get(int(d),0)/64)
    def test_exact_p_superuniform(self):
        for hist in ({1:100},{2:100},{1:12,2:21,3:4}):
            law=m.exact_law(hist)
            for alpha in (.01,.05,.10):self.assertLessEqual(law['pmf'][law['p_exact']<=alpha].sum(),alpha+1e-12)
    def test_checksum_reject(self):
        with tempfile.TemporaryDirectory() as a:
            a=Path(a);fixture(a);(a/'all-guide-counts.tsv').write_text('wrong')
            with self.assertRaises(ValueError):m.load_classes(a)
    def test_repeat_event_information_insufficient(self):
        with tempfile.TemporaryDirectory() as a:
            a=Path(a);fixture(a);obj=json.loads((a/'completion.json').read_text());obj['statistics']['repeated_same_target_events']=1;m.dump(a/'completion.json',obj)
            with self.assertRaises(ValueError):m.load_classes(a)
    def test_bad_counts(self):
        for x in ('-1','1.0','nan',' 1','١',1):
            with self.assertRaises(ValueError):m.natural(x)
    def test_bh(self):
        self.assertEqual(m.bh_count(np.array([[.001,.01,.5],[.1,.2,.3]])).tolist(),[2,0])
    def test_coin_split_variance_monte_carlo(self):
        rng=np.random.default_rng(20260906);D=2*(2*rng.binomial(100,.5,100000)-100)
        self.assertLess(abs(np.var(D)-400),7)
    def test_independent_binomial_law(self):
        law=m.exact_law({1:100});np.testing.assert_allclose(law['pmf'],stats.binom.pmf(np.arange(101),100,.5),atol=1e-15)

if __name__=='__main__':unittest.main(verbosity=2)
