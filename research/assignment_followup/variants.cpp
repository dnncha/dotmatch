// Reuse the validated classification contract; independently preserve raw sequence/quality evidence.
#define main frozen_oracle_main
#include "../assignment_study/oracle.cpp"
#undef main
#include <unordered_set>
struct VariantEvidence { U n=0,qsum=0,q30=0,q40=0,window_qsum=0;int qmin=999,qmax=-1,pos=-1; };
int main(int argc,char**argv){try{
 if(argc!=6)throw runtime_error("usage: variants library.tsv selected.tsv reads.fq.gz output.tsv expected_records");
 Index ix(argv[1]);unordered_map<string,int>byid;for(int i=0;i<(int)ix.g.size();++i)byid[ix.g[i].id]=i;
 ifstream sf(argv[2]);string line;getline(sf,line);unordered_set<int>chosen;
 while(getline(sf,line)){auto row=split(line);if(row.empty()||!byid.count(row[0]))throw runtime_error("Unknown selected guide");chosen.insert(byid[row[0]]);}
 map<pair<int,string>,VariantEvidence>v;string name,s,q;Fastq f(argv[3]);U n=0;
 while(f.next(name,s,q)){n++;auto h=ix.query(s,23);int id=best(h);if(id<0||!chosen.count(id))continue;string w=s.substr(23,L);auto&a=v[{id,w}];a.n++;int diffs=0,pos=-1;for(int i=0;i<L;++i){a.window_qsum+=U(q[23+i]-33);if(w[i]!=ix.g[id].seq[i]){diffs++;pos=i;}}
  if(diffs>1)throw runtime_error("Assignment beyond radius1");a.pos=pos;
  if(pos>=0){int quality=q[23+pos]-33;a.qsum+=quality;a.q30+=quality>=30;a.q40+=quality>=40;a.qmin=min(a.qmin,quality);a.qmax=max(a.qmax,quality);}
 }
 f.close();if(n!=stoull(argv[5]))throw runtime_error("FASTQ record count changed");
 ofstream out(argv[4]);out<<"guide\tgene\treference_sequence\tobserved_sequence\tdistance\tmismatch_position_1based\tcount\tmismatch_phred_sum\tmismatch_q30_reads\tmismatch_q40_reads\tminimum_mismatch_phred\tmaximum_mismatch_phred\twindow_phred_sum\n";
 for(auto&[key,a]:v){int id=key.first;out<<ix.g[id].id<<'\t'<<ix.g[id].gene<<'\t'<<ix.g[id].seq<<'\t'<<key.second<<'\t'<<(a.pos>=0)<<'\t'<<(a.pos+1)<<'\t'<<a.n<<'\t'<<a.qsum<<'\t'<<a.q30<<'\t'<<a.q40<<'\t'<<(a.pos<0?-1:a.qmin)<<'\t'<<a.qmax<<'\t'<<a.window_qsum<<'\n';}
 if(!out)throw runtime_error("Output write failed");cerr<<"Complete "<<n<<" records, "<<chosen.size()<<" selected guides, "<<v.size()<<" observed guide-sequence pairs\n";return 0;
 }catch(const exception&e){cerr<<"ERROR "<<e.what()<<'\n';return 1;}}
