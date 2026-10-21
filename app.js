(function () {
  'use strict';
  const D = window.SCHJERVEN_DATA;
  const $ = id => document.getElementById(id);
  const esc = s => String(s).replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  const fmt = (v, digits=2) => Number.isFinite(v) ? v.toFixed(digits) : 'undefined';
  const signed = v => Number.isFinite(v) ? (v>0?'+':'')+v.toFixed(2) : 'undefined';
  const mean = v => v.reduce((a,b)=>a+b,0)/v.length;
  const ratio = (a,b) => a>0 && b>0 ? Math.log2(a/b) : null;
  const cpms = g => g.counts.map((v,i)=>v*1e6/D.normalization.cpm_denominators[i]);
  const norm = (g, method) => method==='cpm' ? cpms(g) : g.counts.map((v,i)=>v/D.normalization.median_ratio_factors[i]);
  function summarize(g, method) {
    const vals=norm(g,method),cpm=cpms(g);
    return [0,4].map(start=>{
      const control=vals.slice(start,start+2),ik=vals.slice(start+2,start+4);
      const pairs=ik.flatMap(a=>control.map(b=>ratio(a,b)));
      const finite=pairs.filter(v=>v!==null);
      return {effect:ratio(mean(ik),mean(control)),pairs,range:finite.length?[Math.min(...finite),Math.max(...finite)]:null,complete:finite.length===4,
        controlCPM:mean(cpm.slice(start,start+2)),ikCPM:mean(cpm.slice(start+2,start+4))};
    });
  }
  function judge(stats,effect,floor,strict) {
    if(!stats.every(s=>s.controlCPM>=floor&&s.ikCPM>=floor)) return {pass:false,reason:'Below the expression floor'};
    if(!stats.every(s=>s.effect!==null)) return {pass:false,reason:'Undefined mean ratio'};
    const sign=Math.sign(stats[0].effect);
    if(!sign || Math.sign(stats[1].effect)!==sign) return {pass:false,reason:'Different mean directions across models'};
    if(!stats.every(s=>Math.abs(s.effect)>=effect)) return {pass:false,reason:'Mean effect below the chosen floor'};
    if(strict&&!stats.every(s=>s.complete&&s.pairs.every(v=>Math.sign(v)===sign&&Math.abs(v)>=effect))) return {pass:false,reason:'Mean response does not survive every sample comparison'};
    return {pass:true,reason:strict?'Supports a shared readout under the current rule':'Mean-only candidate; inspect the sample comparisons'};
  }
  window.SCHJERVEN_ANALYSIS={cpms,norm,summarize,judge,ratio};
  if(!D || !Array.isArray(D.genes)) { $('gene-title').textContent='Count data could not be loaded. Please reload the page.'; return; }
  const byID=new Map(D.genes.map(g=>[g.id,g])),bySymbol=new Map();
  D.genes.forEach(g=>{if(g.symbol){const k=g.symbol.toUpperCase();if(!bySymbol.has(k))bySymbol.set(k,[]);bySymbol.get(k).push(g);}});
  const cache={cpm:new Map(),median:new Map()};
  const statsFor=(g,m)=>{if(!cache[m].has(g.id))cache[m].set(g.id,summarize(g,m));return cache[m].get(g.id);};
  let selected=(bySymbol.get('ERG')||[])[0]||D.genes[0],queue=[];
  const state=()=>({method:$('normalization').value,effect:Number($('effect-floor').value),floor:Number($('expression-floor').value),strict:$('strict-rule').checked});
  const label=g=>g.symbol||g.id;
  const line=(x1,y1,x2,y2,attr='')=>`<line x1="${x1}" y1="${y1}" x2="${x2}" y2="${y2}" ${attr}/>`;
  const text=(x,y,s,attr='')=>`<text x="${x}" y="${y}" ${attr}>${esc(s)}</text>`;
  function effectsPlot(stats,s) {
    const available=stats.flatMap(t=>[t.effect,...t.pairs]).filter(v=>v!==null);
    const span=Math.max(1,s.effect+0.2,...available.map(Math.abs))*1.13;
    const x=v=>110+(v+span)/(2*span)*460;
    let out='<title id="effect-title">Within-model IK1 versus GFP effects</title><desc id="effect-description">Means and ranges of the four unpaired sample comparisons, in log2 units.</desc>';
    [-span,-span/2,0,span/2,span].forEach(v=>{out+=line(x(v),35,x(v),194,`stroke="${v===0?'#999':'#e6e6e6'}" ${v===0?'stroke-dasharray="3 3"':''}`)+text(x(v),216,fmt(v,1),'text-anchor="middle" font-size="12"');});
    [-s.effect,s.effect].forEach(v=>{out+=line(x(v),35,x(v),194,'stroke="#a7c1c5" stroke-dasharray="2 4"');});
    stats.forEach((t,i)=>{
      const y=75+i*80;out+=text(5,y+4,['MXP5','PDX2'][i],'font-size="15"');
      if(t.range) out+=line(x(t.range[0]),y,x(t.range[1]),y,'stroke="#1f7a8c" stroke-width="3"');
      t.pairs.filter(v=>v!==null).forEach((v,j)=>{out+=`<circle cx="${x(v)}" cy="${y+13+j*2}" r="2.7" fill="#9ebfc5"/>`;});
      if(t.effect!==null) out+=`<circle cx="${x(t.effect)}" cy="${y}" r="6" fill="#1f7a8c"><title>Mean effect ${signed(t.effect)}</title></circle>`;
      out+=text(570,y-17,`mean ${signed(t.effect)}${t.complete?'':' · incomplete range'}`,'text-anchor="end" font-size="12"');
    });
    out+=text(340,247,'Lower ← log₂(IK1 / GFP) → Higher','text-anchor="middle" font-size="12"');
    $('hero-viz').innerHTML=out;
  }
  function samplesPlot(g,s) {
    const cpm=cpms(g),vals=norm(g,s.method),logs=cpm.map(v=>Math.log2(1+v));
    const top=Math.max(1,...logs)*1.1,y=v=>210-v/top*165;let out='<title>Individual library expression</title>';
    for(let i=0;i<=4;i++){const t=top*i/4;out+=line(60,y(t),580,y(t),'stroke="#ddd"')+text(50,y(t)+4,fmt(t,1),'text-anchor="end" font-size="12"');}
    cpm.forEach((v,i)=>{const group=Math.floor(i/2),x=140+group*125+(i%2?-10:10),isIK=group%2===1;out+=`<circle cx="${x}" cy="${y(logs[i])}" r="5" fill="${isIK?'#1f7a8c':'#fff'}" stroke="${isIK?'#1f7a8c':'#555'}" stroke-width="2"><title>${esc(D.samples[i].label||D.samples[i].id)}: ${fmt(v)} CPM</title></circle>`;});
    ['MXP5 GFP','MXP5 IK1','PDX2 GFP','PDX2 IK1'].forEach((v,i)=>out+=text(140+i*125,238,v,'text-anchor="middle" font-size="12"'));
    out+=text(60,24,'log₂(1 + CPM)','font-size="12"')+text(320,268,'Open = GFP · Filled = IK1','text-anchor="middle" font-size="11"');$('sample-viz').innerHTML=out;
    $('sample-table').innerHTML='<table><thead><tr><th>Sample</th><th>Raw count</th><th>CPM</th><th>Effect input</th></tr></thead><tbody>'+D.samples.map((t,i)=>`<tr><td>${esc(t.label||t.id)}</td><td>${g.counts[i]}</td><td>${fmt(cpm[i])}</td><td>${fmt(vals[i])}</td></tr>`).join('')+'</tbody></table>';
  }
  function renderSelected() {
    const s=state(),stats=statsFor(selected,s.method),alt=statsFor(selected,s.method==='cpm'?'median':'cpm');
    const result=judge(stats,s.effect,s.floor,s.strict),other=judge(alt,s.effect,s.floor,s.strict);
    $('gene-title').textContent=`${label(selected)} · ${selected.id}`;effectsPlot(stats,s);samplesPlot(selected,s);
    const sensitivity=result.pass===other.pass?'The shortlist decision is unchanged with the other normalization.':'The shortlist decision changes with normalization; check this before choosing an assay.';
    const rangeText=stats.map((v,i)=>`${['MXP5','PDX2'][i]}: ${v.range?`${signed(v.range[0])} to ${signed(v.range[1])}`:'undefined'}${v.complete?'':' (some comparisons undefined)'}`).join('; ');
    const contextual={HPRT1:'HPRT1 decreases in this dataset: do not assume it is an invariant qPCR reference.',TBP:'A small observed change does not validate TBP as a qPCR reference. Reference stability still needs testing.',MS4A1:'This is an RNA measurement, not evidence of changed surface protein.',IKZF1:'Gene-level IKZF1 counts cannot separate endogenous and induced isoforms.'}[selected.symbol]||'';
    $('readout').innerHTML=`<strong>${esc(result.reason)}.</strong> ${esc(sensitivity)}<br><span>Sample-comparison ranges: ${esc(rangeText)}.</span>${contextual?'<br>'+esc(contextual):''}`;
    const strictNow=judge(stats,s.effect,s.floor,true),strictAlt=judge(alt,s.effect,s.floor,true);
    const next=strictNow.pass&&strictAlt.pass?`For ${label(selected)}, the current evidence supports testing a shared early RNA readout in both models. Repeat it in independent cultures, verify IK1 induction and ask whether the response is retained with an orthogonal measurement.`:`For ${label(selected)}, first resolve the limitation: ${strictNow.pass?'normalization sensitivity':strictNow.reason.toLowerCase()}. Keep the models separate and repeat the RNA measurement with matched processing before treating it as a shared readout.`;
    $('next-measurement').textContent=next+(contextual?' '+contextual:'')+' This is a proposed experiment, not a reported validation.';
    document.querySelectorAll('[data-gene]').forEach(b=>b.setAttribute('aria-pressed',String(b.dataset.gene===selected.id)));
  }
  function renderQueue() {
    const s=state();let meanN=0,strictN=0;queue=[];
    D.genes.forEach(g=>{const stats=statsFor(g,s.method),m=judge(stats,s.effect,s.floor,false),r=judge(stats,s.effect,s.floor,true);if(m.pass)meanN++;if(r.pass)strictN++;if((s.strict?r:m).pass)queue.push({g,stats,rank:s.strict?Math.min(...stats.flatMap(t=>t.pairs.map(Math.abs))):Math.min(...stats.map(t=>Math.abs(t.effect)))});});
    queue.sort((a,b)=>b.rank-a.rank||(a.g.id<b.g.id?-1:a.g.id>b.g.id?1:0));
    $('queue-summary').textContent=`${meanN.toLocaleString()} genes pass the mean-only rule; ${strictN.toLocaleString()} also pass every sample comparison. Showing ${Math.min(queue.length,20)} of ${queue.length.toLocaleString()} current candidates.`;
    $('ranking-note').textContent=s.strict?'Ranked by the weakest of eight sample-comparison effects; largest minimum first. No significance ranking is implied.':'Ranked by the smaller absolute mean effect across the two models. Inspect replicate spread before selecting a readout.';
    const max=queue[0]?.rank||1;
    $('comparison').innerHTML='<div class="comparison-head"><span>Candidate</span><span>Mean log₂ effects · MXP5 / PDX2</span></div>'+(queue.length?queue.slice(0,20).map(q=>`<button type="button" data-gene="${esc(q.g.id)}" aria-pressed="${q.g.id===selected.id}"><span>${esc(label(q.g))}</span><span class="comparison-track"><span class="comparison-fill" style="width:${100*q.rank/max}%"></span></span><span>${signed(q.stats[0].effect)} / ${signed(q.stats[1].effect)}</span></button>`).join(''):'<p class="small">No genes satisfy this combination. An empty shortlist is a result; relax one criterion only if it fits the experimental question.</p>');
    document.querySelectorAll('[data-gene]').forEach(b=>b.addEventListener('click',()=>choose(byID.get(b.dataset.gene))));
    $('comparison-insight').textContent=`At these settings, ${Math.max(0,meanN-strictN).toLocaleString()} mean-only candidates lose support when every observed sample comparison must agree. This is a way to qualify follow-up choices, not to declare those genes biologically unresponsive.`;renderSelected();
  }
  function choose(g){selected=g;$('target-select').value=label(g);$('search-message').textContent='';renderSelected();}
  function correlationTable(){const mat=$('correlation-scale').value==='raw'?D.qc.raw_correlation:D.qc.log_cpm_correlation;const names=D.samples.map((t,i)=>String(i+1));$('correlation-table').innerHTML='<table><caption>Pearson r across all unique gene rows</caption><thead><tr><th>Sample</th>'+names.map(n=>`<th>${n}</th>`).join('')+'</tr></thead><tbody>'+D.samples.map((t,i)=>`<tr><td>${i+1}. ${esc(t.label||t.id)}</td>`+mat[i].map(v=>`<td style="background:rgba(31,122,140,${Math.max(0,v)*0.15})">${fmt(v,3)}</td>`).join('')+'</tr>').join('')+'</tbody></table>';}
  const csv=v=>/[",\n]/.test(String(v))?'"'+String(v).replace(/"/g,'""')+'"':String(v);
  function download(name,content,type){const a=document.createElement('a'),url=URL.createObjectURL(new Blob([content],{type}));a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);}
  $('target-form').addEventListener('submit',e=>{e.preventDefault();const q=$('target-select').value.trim(),id=byID.get(q.toUpperCase()),matches=bySymbol.get(q.toUpperCase())||[];if(id)choose(id);else if(matches.length===1)choose(matches[0]);else if(matches.length>1)$('search-message').textContent='This symbol maps to multiple deposited IDs. Please enter one exact ENSG ID: '+matches.map(g=>g.id).join(', ');else $('search-message').textContent='No exact symbol or gene ID found. Choose a suggestion or search an ENSG ID; no data have been changed.';});
  $('target-select').addEventListener('input',()=>{const q=$('target-select').value.trim().toUpperCase();$('target-options').innerHTML=q.length<2?'':D.genes.filter(g=>g.id.includes(q)||(g.symbol||'').toUpperCase().startsWith(q)).slice(0,15).map(g=>`<option value="${esc(g.symbol||g.id)}">${esc(g.id)}</option>`).join('');});
  ['ERG','MYC','MS4A1','HPRT1','TBP'].filter(k=>bySymbol.has(k)).forEach(k=>{const b=document.createElement('button');b.type='button';b.textContent=k;b.addEventListener('click',()=>choose(bySymbol.get(k)[0]));$('suggestions').append(b);});
  ['normalization','effect-floor','expression-floor','strict-rule'].forEach(id=>$(id).addEventListener('change',renderQueue));$('correlation-scale').addEventListener('change',correlationTable);
  $('download-shortlist').addEventListener('click',()=>{const s=state(),head=['ensembl_id','symbol','normalization','rule','min_abs_log2','min_mean_CPM','MXP5_mean_log2','PDX2_mean_log2','MXP5_range_min','MXP5_range_max','PDX2_range_min','PDX2_range_max','rank_min_abs_effect','MXP5_finite_comparisons','PDX2_finite_comparisons','MXP5_complete_range','PDX2_complete_range'];const rows=queue.map(q=>[q.g.id,q.g.symbol,s.method,s.strict?'every_comparison':'mean_only',s.effect,s.floor,q.stats[0].effect,q.stats[1].effect,...q.stats[0].range,...q.stats[1].range,q.rank,...q.stats.map(t=>t.pairs.filter(v=>v!==null).length),...q.stats.map(t=>t.complete)]);download('schjerven-exploratory-readouts.csv',[head,...rows].map(r=>r.map(csv).join(',')).join('\n')+'\n','text/csv');});
  $('download-manifest').addEventListener('click',()=>{const fields=[...new Set(D.samples.flatMap(s=>Object.keys(s)))];download('schjerven-sample-provenance.csv',[fields,...D.samples.map(s=>fields.map(k=>typeof s[k]==='object'?JSON.stringify(s[k]):s[k]??''))].map(r=>r.map(csv).join(',')).join('\n')+'\n','text/csv');});
  $('download-plan').addEventListener('click',()=>{const s=state();download(`schjerven-${selected.id}-proposed-record.json`,JSON.stringify({status:'PROPOSED_NOT_PERFORMED',source:'GSE310958',gene:{id:selected.id,symbol:selected.symbol},criteria:s,observed_RNA_summary:statsFor(selected,s.method),next_question:$('next-measurement').textContent,controls:['GFP control within each model','Verify IK1 protein induction','Viability and growth state','Validated qPCR reference genes and primer efficiency','Independent biological cultures; technical wells not independent','Matched processing; documented batch; prospective replicate planning'],blank_record:{sample_id:null,model:null,condition:null,biological_replicate:null,batch:null,harvest_time:null,induction:null,processing:null,IK1_protein:null,viability:null,growth:null,qPCR_reference_genes:null,primer_efficiency:null,target_measurement:null,notes:null}},null,2)+'\n','application/json');});
  $('gene-count').textContent=D.genes.length.toLocaleString();renderQueue();correlationTable();
})();
