import test from 'node:test';
import assert from 'node:assert/strict';
import { auditLibrary, parseLibrary, LIMITS } from '../../lib/library-safety.ts';
const panel = (...sequences) => sequences.map((sequence, i) => ({ id: `t${i}`, sequence }));
function run(targets) { const iterator = auditLibrary(targets); for (;;) { const next = iterator.next(); if (next.done) return next.value; } }
function distance(a,b) { let n=0; for(let i=0;i<a.length;i++) if(a[i]!==b[i])n++; return n; }
function mutations(s) { const out=[]; for(let i=0;i<s.length;i++)for(const base of 'ACGT')if(base!==s[i]){ const chars=[...s];chars[i]=base;out.push(chars.join('')); } return out; }
function verify(targets, universe) {
  const report=run(targets), observations=universe ?? [...new Set(targets.flatMap(t=>[t.sequence,...mutations(t.sequence)]))];
  const candidates=s=>targets.filter(t=>distance(t.sequence,s)<=1);
  let distinct=0,ambiguous=0;
  for(const s of observations){const count=candidates(s).length;if(count>0)distinct++;if(count>1)ambiguous++;}
  assert.equal(report.distinct_observations,distinct); assert.equal(report.ambiguous_observations,ambiguous);
  for(const row of report.targets){assert.equal(row.exact_ambiguous_k0,targets.filter(t=>t.sequence===row.sequence).length>1);assert.equal(row.exact_ambiguous_k1,candidates(row.sequence).length>1);assert.equal(row.ambiguous_single_substitutions,mutations(row.sequence).filter(s=>candidates(s).length>1).length);assert.equal(row.possible_single_substitutions,3*row.sequence.length);}
  for(const witness of report.witnesses){assert.equal(witness.candidate_count,candidates(witness.observation).length);assert.deepEqual(witness.target_ids,candidates(witness.observation).map(t=>t.id).sort().slice(0,8));}
  assert.equal(report.targets_with_ambiguous_exact_reads,report.targets.filter(t=>t.exact_ambiguous_k1).length);
  assert.equal(report.targets_with_ambiguous_substitutions,report.targets.filter(t=>t.ambiguous_single_substitutions>0).length);
  return report;
}
test('isolated sequence has 1+3L observations',()=>{const r=verify(panel('AAAAAAAA'));assert.equal(r.distinct_observations,25);assert.equal(r.ambiguous_observations,0);});
test('distance one: four shared observations and ambiguous exact reads',()=>{const r=verify(panel('AAAAAAAA','CAAAAAAA'));assert.equal(r.ambiguous_observations,4);assert.equal(r.targets_with_ambiguous_exact_reads,2);assert.equal(r.targets[0].ambiguous_single_substitutions,3);});
test('distance two: two shared observations, distinct exact reads',()=>{const r=verify(panel('AAAAAAAA','CCAAAAAA'));assert.equal(r.ambiguous_observations,2);assert.equal(r.targets_with_ambiguous_exact_reads,0);});
test('distance three: disjoint radius-one neighbourhoods',()=>{assert.equal(verify(panel('AAAAAAAA','CCCAAAAA')).ambiguous_observations,0);});
test('duplicate sequences remain separate candidates',()=>{const r=verify(panel('AAAAAAAA','AAAAAAAA'));assert.equal(r.ambiguous_observations,25);assert(r.targets.every(t=>t.exact_ambiguous_k0&&t.ambiguous_single_substitutions===24));});
test('three-way overlaps are counted once, not once per pair',()=>verify(panel('AAAAAAAA','CAAAAAAA','GAAAAAAA')));
test('parser supports BOM, CRLF, headers, lowercase, TSV, CSV and blank lines',()=>{assert.deepEqual(parseLibrary('\uFEFFtarget_id\tsequence\r\na\tacgtacgt\r\n\r\nb\ttgcatgca'),[{id:'a',sequence:'ACGTACGT'},{id:'b',sequence:'TGCATGCA'}]);assert.equal(parseLibrary('id,sequence\na,AAAAAAAA')[0].id,'a');assert.equal(parseLibrary('sequence\nacgtacgt')[0].sequence,'ACGTACGT');});
test('malformed input fails closed including an empty ID before a tab',()=>{for(const text of ['', 'id\tsequence', '\tACGTACGT', 'a\tACGTNCGT','a\tACGT','a\tAAAAAAAA\nb\tAAAAAAAAA','a\tAAAAAAAA\na\tCCCCCCCC','a,b,c','"a",AAAAAAAA','a\t'+'A'.repeat(33)])assert.throws(()=>parseLibrary(text),undefined,JSON.stringify(text));});
test('deterministic report is invariant to input row order',()=>{const p=panel('AAAAAAAA','CAAAAAAA','GAAAAAAA');assert.deepEqual(run(p),run([...p].reverse()));});
test('all 65,536 eight-base observations checked against independent distance oracle',()=>{const universe=Array.from({length:4**8},(_,n)=>n.toString(4).padStart(8,'0').replace(/[0-3]/g,d=>'ACGT'[Number(d)]));verify(panel('AAAAAAAA','CAAAAAAA','CCAAAAAA','CCCAAAAA','ACGTACGT','ACGTACGA','AAAAAAAA'),universe);});
test('90 seeded clustered panels at 8,20,32 bases agree with independent oracle',()=>{let state=12345;const rand=()=>{state=(Math.imul(state,1664525)+1013904223)>>>0;return state;};for(const length of [8,20,32])for(let n=0;n<30;n++){const first=Array.from({length},()=> 'ACGT'[rand()%4]).join('');const sequences=[first];for(let j=0;j<7;j++){const chars=[...first];for(let k=0;k<j%4;k++)chars[rand()%length]='ACGT'[rand()%4];sequences.push(chars.join(''));}verify(panel(...sequences));}});
test('maximum supported panel is complete and yields progress',()=>{const p=Array.from({length:LIMITS.targets},(_,i)=>({id:`t${i}`,sequence:i.toString(4).padStart(32,'0').replace(/[0-3]/g,d=>'ACGT'[Number(d)])}));const it=auditLibrary(p);let yields=0;for(;;){const n=it.next();if(n.done){assert.equal(n.value.targets.length,LIMITS.targets);assert(yields>20);break;}assert(n.value.completed<=n.value.total);yields++;}});
test('oversize input is rejected rather than truncated or sampled',()=>{assert.throws(()=>parseLibrary('A'.repeat(LIMITS.characters+1)));assert.throws(()=>run(Array.from({length:2001},(_,i)=>({id:`t${i}`,sequence:'AAAAAAAA'}))));});
test('witness bounds and explicit scientific scope survive serialization',()=>{const r=run(Array.from({length:20},(_,i)=>({id:`t${i}`,sequence:'AAAAAAAA'})));assert.equal(r.witnesses.length,12);assert(r.witnesses_truncated);assert(r.witnesses.every(w=>w.ids_truncated&&w.target_ids.length===8));assert.equal(JSON.parse(JSON.stringify(r)).evidence_level,'exact_combinatorial_audit');assert(r.limitations.some(x=>x.includes('not observed sequencing error rates')));});
