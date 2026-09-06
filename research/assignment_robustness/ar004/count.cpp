// AR004 research counter. Exact candidate enumeration; not a production backend.
#include <cmath>
#include <algorithm>
#include <array>
#include <cctype>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include <zlib.h>
using namespace std;
using U=uint64_t;
static const array<string,5> names={"event_exact","event_best","joint_exact","joint_best","joint_radius"};
vector<string> split(const string&s,char sep){vector<string> v;string x;stringstream z(s);while(getline(z,x,sep))v.push_back(x);return v;}
vector<int> uniq(vector<int> v){sort(v.begin(),v.end());v.erase(unique(v.begin(),v.end()),v.end());return v;}
U mix(U x){x+=0x9e3779b97f4a7c15ULL;x=(x^(x>>30))*0xbf58476d1ce4e5b9ULL;x=(x^(x>>27))*0x94d049bb133111ebULL;return x^(x>>31);}
struct Hit{int pos,id,d;};
struct Call{vector<Hit> hits;array<vector<int>,5> targets;array<vector<int>,3> genes;bool valid=true;};
struct Library{
 vector<string> ids,seq,genes;vector<int> gid;int L=0,cut=0;
 unordered_map<string,vector<int>> left,right;
 Library(const string&path){ifstream f(path);if(!f)throw runtime_error("Cannot open library");string s;getline(f,s);if(s!="id\tsequence\tgene")throw runtime_error("Wrong library header");set<string> seen;map<string,int> gs;
 while(getline(f,s)){if(!s.empty()&&s.back()=='\r')s.pop_back();auto r=split(s,'\t');if(r.size()!=3||r[0].empty()||r[2].empty()||!seen.insert(r[0]).second)throw runtime_error("Invalid or duplicate reference identity");if(!L){L=r[1].size();cut=L/2;}if(L<2||r[1].size()!=size_t(L)||r[1].find_first_not_of("ACGTRYSWKMBDHVN")!=string::npos)throw runtime_error("Invalid reference sequence");int i=seq.size();ids.push_back(r[0]);seq.push_back(r[1]);auto [it,added]=gs.emplace(r[2],gs.size());if(added)genes.push_back(r[2]);gid.push_back(it->second);left[r[1].substr(0,cut)].push_back(i);right[r[1].substr(cut)].push_back(i);}
 if(seq.empty())throw runtime_error("Empty library");}
 vector<pair<int,int>> candidates(const string&w)const{vector<int> ids;auto a=left.find(w.substr(0,cut));if(a!=left.end())ids=a->second;auto b=right.find(w.substr(cut));if(b!=right.end())ids.insert(ids.end(),b->second.begin(),b->second.end());ids=uniq(ids);vector<pair<int,int>> out;for(int i:ids){int d=0;for(int j=0;j<L&&d<=1;j++)d+=w[j]!=seq[i][j];if(d<=1)out.emplace_back(i,d);}return out;}
 Call decode(const string&s,const vector<int>&positions)const{Call c;if(s.size()<size_t(positions.back()+L)||s.find_first_not_of("ACGTRYSWKMBDHVN")!=string::npos){c.valid=false;return c;}
 for(int p:positions){auto hits=candidates(s.substr(p,L));vector<int> exact,best;int bd=2;for(auto [i,d]:hits){c.hits.push_back({p,i,d});if(!d)exact.push_back(i);bd=min(bd,d);}for(auto [i,d]:hits)if(d==bd)best.push_back(i);if(exact.size()==1)c.targets[0].push_back(exact[0]);if(best.size()==1)c.targets[1].push_back(best[0]);}
 int bd=2;for(auto h:c.hits)bd=min(bd,h.d);for(auto h:c.hits){if(!h.d)c.targets[2].push_back(h.id);if(h.d==bd)c.targets[3].push_back(h.id);c.targets[4].push_back(h.id);}for(int k=2;k<5;k++){c.targets[k]=uniq(c.targets[k]);for(int i:c.targets[k])c.genes[k-2].push_back(gid[i]);c.genes[k-2]=uniq(c.genes[k-2]);}return c;}
};
struct Reader{gzFile f;Reader(const string&p){f=gzopen(p.c_str(),"rb");if(!f)throw runtime_error("Cannot open FASTQ");gzbuffer(f,1<<20);}~Reader(){if(f)gzclose(f);}bool line(string&s){s.clear();char b[8192];while(true){char*p=gzgets(f,b,sizeof b);if(!p){int err;gzerror(f,&err);if(err!=Z_OK&&err!=Z_STREAM_END)throw runtime_error("Gzip input error");return !s.empty();}s+=b;if(s.size()>1000000)throw runtime_error("Oversize FASTQ line");if(!s.empty()&&s.back()=='\n'){s.pop_back();if(!s.empty()&&s.back()=='\r')s.pop_back();return true;}}}bool next(string&seq){string h,p,q;if(!line(h))return false;if(!line(seq)||!line(p)||!line(q)||h.size()<2||h[0]!='@'||p.empty()||p[0]!='+'||seq.empty()||seq.size()!=q.size())throw runtime_error("Malformed FASTQ");for(unsigned char c:q)if(c<33||c>126)throw runtime_error("Invalid quality");for(char&c:seq)c=toupper((unsigned char)c);return true;}};
int main(int argc,char**argv){try{
 if(argc<4){throw runtime_error("usage: count LIB OFFSETS --probe | FASTQ OUT RATE SEED CACHE_LIMIT");}
 Library lib(argv[1]);vector<int> pos;for(auto&s:split(argv[2],',')){size_t n=0;int p=stoi(s,&n);if(n!=s.size()||p<0)throw runtime_error("Invalid offset");pos.push_back(p);}sort(pos.begin(),pos.end());if(pos.empty()||adjacent_find(pos.begin(),pos.end())!=pos.end())throw runtime_error("Invalid positions");
 if(string(argv[3])=="--probe"){string s;while(getline(cin,s)){auto c=lib.decode(s,pos);cout<<(c.valid?1:0)<<'\t';for(auto h:c.hits)cout<<h.pos<<','<<h.id<<','<<h.d<<';';for(auto&t:c.targets){cout<<'\t';for(int i:t)cout<<i<<',';}cout<<'\n';}return 0;}
 if(argc<8){throw runtime_error("Missing arguments");}
 filesystem::path out=argv[4];if(filesystem::exists(out))throw runtime_error("Refuse existing output");filesystem::create_directories(out);
 double rate=stod(argv[5]);if(!isfinite(rate)||rate<0||rate>1)throw runtime_error("Invalid thinning rate");U seed=stoull(argv[6]);size_t cap=stoull(argv[7]);
 array<vector<U>,10> counts;for(auto&v:counts)v.assign(lib.seq.size(),0);array<vector<U>,3> lower,upper;for(int i=0;i<3;i++){lower[i].assign(lib.genes.size(),0);upper[i].assign(lib.genes.size(),0);}array<array<U,4>,5> qc{};array<array<U,4>,3> gqc{};array<U,2> events{},matched{},extra{},multireads{},same_gene{},cross_gene{};
 map<pair<int,int>,U> pairs;unordered_map<string,Call> cache;cache.reserve(cap?cap:1);U N=0,thinN=0;Reader f(argv[3]);string seq;
 while(f.next(seq)){N++;bool thin=rate>=1||(double(mix(seed^N)>>11)*(1.0/9007199254740992.0)<rate);thinN+=thin;
 string key=seq;
 Call value;const Call*c;auto found=cache.find(key);if(found!=cache.end())c=&found->second;else{value=lib.decode(seq,pos);if(cap){if(cache.size()>=cap)cache.clear();auto it=cache.emplace(move(key),move(value));c=&it.first->second;}else c=&value;}
 for(int k=0;k<5;k++){auto&t=c->targets[k];int st=!c->valid?3:t.empty()?2:((k<2||t.size()==1)?0:1);qc[k][st]++;if(c->valid&&(k<2||t.size()==1))for(int i:t){counts[k][i]++;if(thin)counts[k+5][i]++;}}
 for(int k=0;k<3;k++){auto&g=c->genes[k];int st=!c->valid?3:g.empty()?2:g.size()==1?0:1;gqc[k][st]++;if(c->valid){if(g.size()==1)lower[k][g[0]]++;for(int i:g)upper[k][i]++;}}
 for(int k=0;k<2;k++){auto&t=c->targets[k];if(!t.empty()){matched[k]++;events[k]+=t.size();extra[k]+=t.size()-1;if(t.size()>1){multireads[k]++;vector<int> genes;for(int i:t)genes.push_back(lib.gid[i]);if(uniq(genes).size()==1)same_gene[k]++;else cross_gene[k]++;}}}
 auto ids=uniq(c->targets[1]);for(size_t i=0;i<ids.size();i++)for(size_t j=i+1;j<ids.size();j++)pairs[{ids[i],ids[j]}]++;
 if(N%10000000==0)cerr<<N<<" records\n";
 }
 ofstream f1(out/"guide-counts.tsv");f1<<"id\tgene";for(auto&s:names)f1<<'\t'<<s;for(auto&s:names)f1<<"\tthin_"<<s;f1<<'\n';for(size_t i=0;i<lib.seq.size();i++){f1<<lib.ids[i]<<'\t'<<lib.genes[lib.gid[i]];for(auto&v:counts)f1<<'\t'<<v[i];f1<<'\n';}
 ofstream fg(out/"gene-counts.tsv");fg<<"gene";for(int k=2;k<5;k++)fg<<'\t'<<names[k]<<"_lower\t"<<names[k]<<"_upper";fg<<'\n';for(size_t g=0;g<lib.genes.size();g++){fg<<lib.genes[g];for(int k=0;k<3;k++)fg<<'\t'<<lower[k][g]<<'\t'<<upper[k][g];fg<<'\n';}
 ofstream fq(out/"qc.tsv");fq<<"policy\tresolution\tunique\tambiguous\tnone\tinvalid\n";for(int k=0;k<5;k++){fq<<names[k]<<(k<2?"\tevent_matched_record":"\tguide");for(U v:qc[k])fq<<'\t'<<v;fq<<'\n';}for(int k=0;k<3;k++){fq<<names[k+2]<<"\tgene";for(U v:gqc[k])fq<<'\t'<<v;fq<<'\n';}
 ofstream fp(out/"shared-read-pairs.tsv");fp<<"guide_a\tguide_b\tgene_a\tgene_b\tshared_records\n";for(auto&[p,n]:pairs)fp<<lib.ids[p.first]<<'\t'<<lib.ids[p.second]<<'\t'<<lib.genes[lib.gid[p.first]]<<'\t'<<lib.genes[lib.gid[p.second]]<<'\t'<<n<<'\n';
 ofstream fm(out/"metrics.json");fm<<"{\"records\":"<<N<<",\"thinned_records\":"<<thinN<<",\"policies\":{";for(int k=0;k<2;k++){if(k)fm<<',';if(events[k]!=matched[k]+extra[k])throw runtime_error("Count conservation failed");fm<<'"'<<names[k]<<"\":{\"count_events\":"<<events[k]<<",\"matched_records\":"<<matched[k]<<",\"extra_events\":"<<extra[k]<<",\"multiply_counted_records\":"<<multireads[k]<<",\"within_gene_multiple\":"<<same_gene[k]<<",\"cross_gene_multiple\":"<<cross_gene[k]<<'}';}fm<<"}}\n";
 for(int k=0;k<5;k++){U sum=0;for(U v:qc[k])sum+=v;if(sum!=N)throw runtime_error("Read budget failed");}
 f1.close();fg.close();fq.close();fp.close();fm.close();
 if(!f1||!fg||!fq||!fp||!fm){throw runtime_error("Output write failure");}
 return 0;
 }catch(const exception&e){cerr<<"ERROR: "<<e.what()<<'\n';return 1;}}
