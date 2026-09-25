// Independent exhaustive-neighborhood Hamming oracle and per-read accounting.
// No DotMatch or guide-counter library/API is used here.
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <tuple>
#include <vector>
#include <zlib.h>
using namespace std;
using U=uint64_t;
constexpr int L=19;
struct Guide { string id,seq,gene; U code; };
struct Hit { int exact=-1, first=-1, n=0; };
int base(char c) {switch(c){case 'A':return 0;case 'C':return 1;case 'G':return 2;case 'T':return 3;default:return -1;}}
U encode(const string&s){U v=0;for(char c:s){int b=base(c);if(b<0)throw runtime_error("Non-ACGT library");v=(v<<2)|b;}return v;}
vector<string> split(const string&s,char sep='\t'){vector<string>r;string x;istringstream ss(s);while(getline(ss,x,sep))r.push_back(x);return r;}
struct Fastq {
 gzFile f; string path; U ordinal=0; char buf[8192];
 Fastq(const string&p):path(p){f=gzopen(p.c_str(),"rb");if(!f)throw runtime_error("Cannot open "+p);gzbuffer(f,1<<20);}
 bool line(string&s){s.clear();while(true){char*p=gzgets(f,buf,sizeof(buf));if(!p){int e;string m=gzerror(f,&e);if(e!=Z_OK&&e!=Z_STREAM_END)throw runtime_error(m);return !s.empty();}s+=p;if(!s.empty()&&s.back()=='\n'){s.pop_back();if(!s.empty()&&s.back()=='\r')s.pop_back();return true;}}}
 bool next(string&name,string&seq,string&q){string plus;if(!line(name))return false;if(!line(seq)||!line(plus)||!line(q)||name.empty()||name[0]!='@'||plus.empty()||plus[0]!='+'||seq.size()!=q.size())throw runtime_error("Malformed FASTQ "+path+":"+to_string(ordinal+1));for(char c:q)if(c<33||c>126)throw runtime_error("Invalid quality");++ordinal;return true;}
 void close(){if(f){int e=gzclose(f);f=nullptr;if(e!=Z_OK)throw runtime_error("Gzip close failure");}}
 ~Fastq(){if(f)gzclose(f);}
};
struct Index {
 vector<Guide>g;unordered_map<U,Hit>m;unordered_map<U,int> exact;
 Index(string file){ifstream f(file);if(!f)throw runtime_error("library missing");string l;getline(f,l);set<string>ids;while(getline(f,l)){if(l.empty())continue;auto a=split(l);if(a.size()!=3||a[1].size()!=L||!ids.insert(a[0]).second)throw runtime_error("Invalid library");Guide t{a[0],a[1],a[2],encode(a[1])};if(exact.count(t.code))throw runtime_error("duplicate sequence");exact[t.code]=g.size();g.push_back(t);}if(g.empty())throw runtime_error("empty library");m.reserve(g.size()*58);for(int i=0;i<(int)g.size();++i){U v=g[i].code;auto& h=m[v];h.exact=i;if(!h.n)h.first=i;h.n++;for(int p=0;p<L;++p){int sh=2*p,b=(v>>sh)&3;for(int x=0;x<4;++x)if(x!=b){U w=(v&~(U(3)<<sh))|(U(x)<<sh);auto&j=m[w];if(!j.n)j.first=i;j.n++;}}}cerr<<"Oracle indexed "<<g.size()<<" guides / "<<m.size()<<" words\n";}
 Hit query(const string&s,int start,bool allow_non_acgt=true)const{if(start<0||start+L>(int)s.size())return {};U v=0;int bad=0,bp=-1;for(int p=0;p<L;++p){int b=base(s[start+p]);v<<=2;if(b<0){bad++;bp=L-1-p;}else v|=b;}if(!bad){auto it=m.find(v);return it==m.end()?Hit{}:it->second;}if(!allow_non_acgt||bad>1)return {};Hit h;for(int b=0;b<4;++b){auto it=exact.find(v|(U(b)<<(2*bp)));if(it!=exact.end()){if(!h.n)h.first=it->second;h.n++;}}return h;}
 Hit brute(const string&s,int start)const{Hit h;if(start<0||start+L>(int)s.size())return h;for(int i=0;i<(int)g.size();++i){int d=0;for(int p=0;p<L;++p){d+=s[start+p]!=g[i].seq[p];if(d>1)break;}if(d<=1){if(!h.n)h.first=i;h.n++;if(!d)h.exact=i;}}return h;}
};
int best(Hit h){return h.exact>=0?h.exact:(h.n==1?h.first:-1);}
int radius(Hit h){return h.n==1?h.first:-1;}
int state(Hit h,int p,bool invalid=false){if(invalid)return 3;if(p==0)return h.exact>=0?0:2;if(p==1)return h.n==1?0:h.n?1:2;return best(h)>=0?0:h.n?1:2;}
string joinints(const vector<int>&v){string r;for(int x:v){if(!r.empty())r+=",";r+=to_string(x);}return r;}
int main(int argc,char**argv){try{
 if(argc!=5)throw runtime_error("usage: oracle library.tsv reads.fastq.gz prefix calibration_records");
 Index ix(argv[1]);string path=argv[2],pre=argv[3];U calibration=stoull(argv[4]);
 array<vector<U>,2>hist;for(auto&h:hist)h.resize(1000);string name,s,q;U caln=0;
 {Fastq f(path);while(caln<calibration&&f.next(name,s,q)){if(s.size()>1000)throw runtime_error("long read outside study contract");for(int o=0;o+L<=(int)s.size();++o){auto h=ix.query(s,o,false);if(h.exact>=0)hist[0][o]++;if(best(h)>=0)hist[1][o]++;}caln++;}}
 array<vector<int>,2>offsets;ofstream of(pre+".offsets.tsv");of<<"mode\toffset\tcalibration_hits\tcalibration_total_window_hits\tretained\n";
 for(int p=0;p<2;++p){U total=0;for(U x:hist[p])total+=x;if(!total)throw runtime_error("No calibration matches");for(int o=0;o<1000;++o)if(hist[p][o]){bool keep=double(hist[p][o])/double(total)>=0.0025;if(keep)offsets[p].push_back(o);of<<(p?"best":"exact")<<'\t'<<o<<'\t'<<hist[p][o]<<'\t'<<total<<'\t'<<keep<<'\n';}}
 vector<string>pol={"fixed_exact","fixed_radius1","fixed_best1","multi_exact","distinct_exact","multi_best1","distinct_best1","joint_best1"};
 vector<vector<U>>counts(pol.size(),vector<U>(ix.g.size())),held=counts;
 array<array<U,4>,3>states{};array<map<string,U>,2>qc;
 array<map<string,U>,3>trans;map<tuple<int,int,int>,U> edges;
 gzFile ledger=gzopen((pre+".multihit.tsv.gz").c_str(),"wb1");if(!ledger)throw runtime_error("ledger open");gzputs(ledger,"record_ordinal\tread_id\tmode\toffsets\tguide_ids\texcess_same_guide_contributions\tdistinct_guides\n");
 ofstream examples(pre+".examples.tsv");examples<<"record_ordinal\tsequence\tmode\toffsets\tguide_ids\n";
 ofstream checks(pre+".exhaustive_checks.tsv");checks<<"record_ordinal\toffset\tn_radius\texact_id\tbest_id\n";
 array<int,4>strata{};U checked=0,n=0,nonacgt=0;int examplecount=0;
 Fastq f(path);
 while(f.next(name,s,q)){
  ++n;bool invalid=s.size()<23+L;auto fh=ix.query(s,23);array<int,3> fid={fh.exact,radius(fh),best(fh)};array<int,3>st;
  for(int p=0;p<3;++p){st[p]=state(fh,p,invalid);states[p][st[p]]++;if(fid[p]>=0){counts[p][fid[p]]++;if(n>caln)held[p][fid[p]]++;}}
  for(int p=0;p<3;++p){int a=p==2?1:0,b=p==0?1:2;trans[p][to_string(st[a])+"_"+to_string(st[b])]++;}
  bool bad=false;if(!invalid)for(int j=23;j<23+L;++j)bad|=base(s[j])<0;nonacgt+=bad;
  int str=fh.exact>=0?0:fh.n==1?1:fh.n>1?2:3;
  bool test=(!invalid&&((n>caln&&n%10007==0&&checked<1000)||(n>caln&&strata[str]<20)));
  if(test){auto b=ix.brute(s,23);if(b.n!=fh.n||b.exact!=fh.exact||best(b)!=best(fh))throw runtime_error("Exhaustive fixed oracle discrepancy");strata[str]++;checked++;checks<<n<<"\t23\t"<<b.n<<'\t'<<(b.exact<0?".":ix.g[b.exact].id)<<'\t'<<(best(b)<0?".":ix.g[best(b)].id)<<'\n';}
  for(int mode=0;mode<2;++mode){
   vector<int>ids,os;set<int>unique;bool jointamb=false;int jointdist=2,jointid=-1;
   for(int o:offsets[mode]){
    auto h=ix.query(s,o,false);int id=mode?best(h):h.exact;
    if(mode&&id>=0&&o!=23&&best(fh)>=0)edges[{best(fh),id,o}]++;
    if(id>=0){ids.push_back(id);os.push_back(o);unique.insert(id);int p=mode?5:3;counts[p][id]++;if(n>caln)held[p][id]++;}
    if(mode&&h.n){int dist=h.exact>=0?0:1;int candidate=h.exact>=0?h.exact:h.first;bool amb=h.exact<0&&h.n>1;if(dist<jointdist){jointdist=dist;jointid=candidate;jointamb=amb;}else if(dist==jointdist){if(amb||candidate!=jointid)jointamb=true;}}
   }
   int p=mode?6:4;for(int id:unique){counts[p][id]++;if(n>caln)held[p][id]++;}
   if(mode&&jointdist<=1&&!jointamb){counts[7][jointid]++;if(n>caln)held[7][jointid]++;}
   set<string> genes;for(int id:unique)genes.insert(ix.g[id].gene);
   auto&c=qc[mode];c["same_gene_cross_guide_reads"]+=(unique.size()>1&&genes.size()==1);c["cross_gene_reads"]+=genes.size()>1;c["reads"]++;c["matched_reads"]+=!ids.empty();c["window_contributions"]+=ids.size();c["distinct_guide_contributions"]+=unique.size();c["same_guide_excess"]+=ids.size()-unique.size();c["multiple_hit_reads"]+=ids.size()>1;c["cross_guide_reads"]+=unique.size()>1;c["same_only_multiple_reads"]+=(ids.size()>1&&unique.size()==1);
   if(n>caln){c["heldout_reads"]++;c["heldout_window_contributions"]+=ids.size();c["heldout_same_guide_excess"]+=ids.size()-unique.size();}
   if(ids.size()>1){string gs;for(int id:ids){if(!gs.empty())gs+=",";gs+=ix.g[id].id;}string line=to_string(n)+"\t"+name.substr(1)+"\t"+(mode?"best":"exact")+"\t"+joinints(os)+"\t"+gs+"\t"+to_string(ids.size()-unique.size())+"\t"+to_string(unique.size())+"\n";if(gzputs(ledger,line.c_str())<0)throw runtime_error("ledger write");if(examplecount<100){examples<<n<<'\t'<<s<<'\t'<<(mode?"best":"exact")<<'\t'<<joinints(os)<<'\t'<<gs<<'\n';examplecount++;}}
  }
  if(n%5000000==0)cerr<<pre<<" processed "<<n<<"\n";
 }
 f.close();if(gzclose(ledger)!=Z_OK)throw runtime_error("ledger close");
 for(int hold=0;hold<2;++hold){ofstream out(pre+(hold?".heldout.counts.tsv":".counts.tsv"));out<<"guide\tgene";for(auto p:pol)out<<'\t'<<p;out<<'\n';auto&cs=hold?held:counts;for(int i=0;i<(int)ix.g.size();++i){out<<ix.g[i].id<<'\t'<<ix.g[i].gene;for(auto&v:cs)out<<'\t'<<v[i];out<<'\n';}}
 ofstream edgeout(pre+".offset_edges.tsv");edgeout<<"canonical_guide\talternate_guide\toffset\tcanonical_gene\talternate_gene\tread_count\n";for(auto [key,v]:edges){auto [a,b,o]=key;edgeout<<ix.g[a].id<<'\t'<<ix.g[b].id<<'\t'<<o<<'\t'<<ix.g[a].gene<<'\t'<<ix.g[b].gene<<'\t'<<v<<'\n';}
 ofstream stats(pre+".qc.tsv");stats<<"mode\tmetric\tvalue\n";stats<<"all\treads\t"<<n<<"\nall\tcalibration_records\t"<<caln<<"\nall\texhaustive_checks\t"<<checked<<"\nall\tfixed_non_acgt\t"<<nonacgt<<'\n';
 for(int p=0;p<3;++p)for(int t=0;t<4;++t)stats<<pol[p]<<'\t'<<vector<string>{"unique","ambiguous","unmatched","invalid"}[t]<<'\t'<<states[p][t]<<'\n';
 for(int p=0;p<2;++p)for(auto[k,v]:qc[p])stats<<(p?"multi_best1":"multi_exact")<<'\t'<<k<<'\t'<<v<<'\n';
 for(int p=0;p<3;++p)for(auto[k,v]:trans[p])stats<<"transition_"<<p<<'\t'<<k<<'\t'<<v<<'\n';
 cerr<<"COMPLETE "<<pre<<" "<<n<<" records, "<<checked<<" exhaustive checks\n";return 0;
 }catch(const exception&e){cerr<<"ERROR "<<e.what()<<"\n";return 1;}}
