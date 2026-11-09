import fs from 'node:fs/promises';
import path from 'node:path';
import {Workbook, SpreadsheetFile} from '@oai/artifact-tool';

const base=path.dirname(new URL(import.meta.url).pathname);
const payload=JSON.parse(await fs.readFile(path.join(base,'workbook_payload.json'),'utf8'));
const previewOnly=process.argv.includes('--preview-only');
if(previewOnly)payload.tables=payload.tables.filter(t=>!t.name.includes('Risk'));
const output='/Users/4482173/Documents/Projects/IMO_workshop_14/tMN/outputs/01a1ac40-905d-7660-b1ab-5500bfe0661f/tmn_all_genes_gt2_results.xlsx';
const wb=Workbook.create();
const navy='#203B59', pale='#EFF3F7', ink='#172B40';
const label={patient_id:'Patient ID',gene:'Gene',s2:'Fitted s₂',median_s2:'Median s₂',n_fitted_patient_gene_groups:'Fitted patients',median_cluster_bootstrap_95_low:'Median bootstrap 95% lower',median_cluster_bootstrap_95_high:'Median bootstrap 95% upper',n_zero_boundary:'s₂ = 0 boundary',n_mismatch_flag:'Fit mismatch flag',exploratory_small_n:'Exploratory n < 5',mu:'μ per division',r_per_day:'Additional growth scale per day',D0:'Baseline double-hit cells',day:'Day from CAR T',P_first_passage:'First-passage probability (0–1)',first_vaf_percent:'First VAF (%)',last_vaf_percent:'Last VAF (%)',predicted_last_vaf_percent:'Predicted last VAF (%)',tmn_source_flag:'tMN source flag (0/1)',source_excel_row:'Endpoint source row',gene_retained:'Gene retained (>2 eligible patients)',eligible_endpoint_pair:'Eligible endpoint pair',analysis_included:'Endpoint included',retained_by_gt2_filter:'Retained (>2 eligible patients)',eligible_patient_gene_groups:'Eligible patients',all_patient_gene_groups:'All patients',eligible_variant_records:'Eligible variants',all_variant_records:'All variants',first_vaf_source_value:'First VAF original source',last_vaf_source_value:'Last VAF original source',tmn_source_excel_rows:'tMN-positive source rows',median_difference_s2:'Difference in median s₂',bonferroni_familywise_95_low:'Bonferroni familywise 95% lower',bonferroni_familywise_95_high:'Bonferroni familywise 95% upper'};
function col(n){let s='';for(n++;n>0;n=Math.floor((n-1)/26))s=String.fromCharCode(65+(n-1)%26)+s;return s;}
function formatFor(k){
 if(k==='P_first_passage')return '0.000E+00';
 if(k.includes('D_fraction')||k==='meanfield_endpoint_D_fraction')return '0.000E+00';
 if(k==='mu'||k.includes('mass_error')||k.includes('solver_'))return '0.00E+00';
 if(k==='founder_count')return '#,##0.00';
 if(k==='output_roundoff_bound_adjustments')return '#,##0';
 if(/(^n_|draws|variants|variant_records|patient_gene_groups|patients|count|cells|source_excel_row|^D0$|^day$|_day$|^tmn_source_flag$|planned_contrasts)/.test(k)&&!k.includes('r_per_day')&&!k.includes('growth'))return '#,##0';
 return '0.0000';
}
function width(k){
 if(k==='gene')return 13;
 if(k==='patient_id')return 17;
 if(k==='variant_id')return 31;
 if(k==='scenario')return 38;
 if(k.includes('flags')||k.includes('caution'))return 75;
 if(k.includes('interpretation')||k==='status'||k.includes('bootstrap_status'))return 58;
 if(k==='contrast')return 32;
 if(k==='parameter')return 55;
 if(k==='scenario_role')return 46;
 if(k==='units')return 28;
 return 22;
}
const manifest=[];
for(const t of payload.tables){
 const s=wb.worksheets.add(t.name); const start=t.name==='Summary'?6:1;
 const matrix=[t.columns.map(k=>label[k]??k.replaceAll('_',' ')),...t.rows];
 const last=start+matrix.length-1;
 s.getRangeByIndexes(start-1,0,matrix.length,t.columns.length).values=matrix;
 s.showGridLines=false;
 s.getRangeByIndexes(start-1,0,matrix.length,t.columns.length).format.font={name:'Arial',size:10,color:ink};
 s.getRangeByIndexes(start-1,0,matrix.length,t.columns.length).format.verticalAlignment='center';
 s.getRangeByIndexes(start,0,t.rows.length,t.columns.length).format.rowHeight=20;
 for(let j=0;j<t.columns.length;j++){
  const k=t.columns[j]; const range=s.getRange(`${col(j)}${start}:${col(j)}${last}`);
  range.format.columnWidth=width(k);
  if(t.kinds[k]==='number'){
   s.getRange(`${col(j)}${start+1}:${col(j)}${last}`).setNumberFormat(formatFor(k));
   s.getRange(`${col(j)}${start+1}:${col(j)}${last}`).format.horizontalAlignment='right';
  }
  if(t.kinds[k]==='text')s.getRange(`${col(j)}${start+1}:${col(j)}${last}`).setNumberFormat('@');
  if(k.includes('flags')||k.includes('caution'))s.getRange(`${col(j)}${start+1}:${col(j)}${last}`).format.wrapText=true;
 }
 const longColumns=t.columns.map((k,j)=>(k.includes('flags')||k.includes('caution'))?j:-1).filter(j=>j>=0);
 if(longColumns.length){
  const heights=t.rows.map(row=>Math.max(1,...longColumns.map(j=>Math.ceil(String(row[j]??'').length/Math.floor(width(t.columns[j])*.95)))));
  for(let i=0;i<heights.length;){
   let end=i+1;while(end<heights.length&&heights[end]===heights[i])end++;
   if(heights[i]>1)s.getRangeByIndexes(start+i,0,end-i,t.columns.length).format.rowHeight=16*heights[i]+4;
   i=end;
  }
 }
 const h=s.getRangeByIndexes(start-1,0,1,t.columns.length);
 h.format={fill:navy,font:{name:'Arial',size:10,bold:true,color:'#FFFFFF'},horizontalAlignment:'center',verticalAlignment:'center',wrapText:true,rowHeight:42};
 const table=s.tables.add(`A${start}:${col(t.columns.length-1)}${last}`,true,t.name.replaceAll(' ','')+'Table');
 table.showFilterButton=true;
 table.style='TableStyleMedium2';
 if(t.name!=='Summary'){s.freezePanes.freezeRows(start);s.freezePanes.freezeColumns(t.columns[0]==='patient_id'||t.columns[0]==='scenario'?2:1);}
 if(t.name==='Patient Gene Fits'){
   const j=t.columns.indexOf('boundary');
   s.getRange(`${col(j)}2:${col(j)}${last}`).conditionalFormats.add('cellIs',{operator:'equal',formula:'TRUE',format:{fill:'#FFF0D0',font:{color:'#7B4F00'}}});
 }
 manifest.push({name:t.name,source:t.source,start,columns:t.columns,rows:t.rows.length});
 console.log(`Populated ${t.name}: ${t.rows.length} rows`);
}
const summary=wb.worksheets.getItem('Summary'); summary.tabColor=navy;
summary.getRange('A2').values=[['Moran model: six-gene endpoint study']];
summary.getRange('A2').format.font={name:'Arial',size:14,bold:true,color:navy};
summary.getRange('A2:H2').format.rowHeight=25;
summary.getRange('A3').values=[['Genes retained with >2 eligible patients. Reference kinetics fixed; patient/gene s₂ constrained to ≥0.']];
summary.getRange('A4').values=[['Imported scientific results. Refresh with the analysis scripts; Excel does not refit or simulate the model.']];
summary.getRange('A3:H4').format.font={name:'Arial',size:10,italic:true,color:'#46586B'};
summary.getRange('A14').values=[['Data coverage']]; summary.getRange('A14').format.font={name:'Arial',size:11,bold:true,color:navy};
const stats=[['All endpoint variant records',"=COUNTA('Endpoints'!A2:A172)"],['Selected patient/gene groups',"=COUNTA('Patient Gene Fits'!A2:A109)"],['Fitted patient/gene groups',"=COUNTIFS('Patient Gene Fits'!F2:F109,\"conditional_descriptive_fit\")"],['Eligible selected variant records',"=SUM('Patient Gene Fits'!D2:D109)"],['s₂ = 0 boundary fits','=SUM(F7:F12)'],['Fits with mismatch flag','=SUM(G7:G12)']];
for(let i=0;i<stats.length;i++){summary.getRange(`A${15+i}`).values=[[stats[i][0]]];summary.getRange(`D${15+i}`).formulas=[[stats[i][1]]];}
summary.getRange('A15:D20').format.font={name:'Arial',size:10,color:ink};summary.getRange('D15:D20').setNumberFormat('#,##0');
const summaryNotes=[
 'Bootstrap intervals describe patient resampling conditional on the reference model and endpoint selection.',
 'ASXL1 (n=3) and CBL (n=4) comparisons are exploratory. Boundary estimates are retained in summaries.',
 'The 15 planned gene contrasts include nominal and Bonferroni-adjusted percentile intervals.',
 'Risk tables report conditional first passage to ≥1,001 aggregate double-hit cells, using fixed fitted s₂.',
 'Risk scenarios are exploratory model calculations; their envelopes are not confidence intervals.',
 'Source flag 0 means no tMN record in the designated sheet. It is not a verified clinical negative.'
];
summaryNotes.forEach((v,i)=>summary.getRange(`A${23+i}`).values=[[v]]);
summary.getRange('A23:H28').format.font={name:'Arial',size:10,color:ink};
const endpoints=wb.worksheets.getItem('Endpoints');
endpoints.getRange('S1:T4').values=[['Source','Citation'],['Endpoint observations','41588_2026_2795_MOESM9_ESM.xlsx; Gene_Level_Expansion_Data; source_excel_row'],['tMN status','41588_2026_2795_MOESM7_ESM.xlsx; tMN_VAF; Paper_ID joined to Patient_ID; actual header Timpoint = tMN'],['Original missing value','Patient_103 / ASXL1 baseline is source text NA. Its numeric VAF remains blank.']];
endpoints.getRange('S1:S4').format.columnWidth=25;endpoints.getRange('T1:T4').format.columnWidth=115;
endpoints.getRange('S1:T4').format.font={name:'Arial',size:10,color:ink};
const outcome=wb.worksheets.getItem('Outcome Evidence');
outcome.getRange('K1:L4').values=[['Source','Interpretation'],['Clinical timing','41588_2026_2795_MOESM7_ESM.xlsx; Auto_CAR_tMN; clinical_source_excel_row'],['Source-status join','Patient_ID matched to Paper_ID in tMN_VAF; no source record does not establish a clinical negative'],['Time origins','Patients 18 and 69 have multiple CAR timing rows. Timing origins require reconciliation before clinical risk calibration.']];
outcome.getRange('K1:K4').format.columnWidth=25;outcome.getRange('L1:L4').format.columnWidth=115;
outcome.getRange('K1:L4').format.font={name:'Arial',size:10,color:ink};
const settings=wb.worksheets.getItem('Model Settings'); settings.tabColor='#6483A1';
const notes=[
 ['Model','Fixed-N Moran birth/death replacement with event clock N × b. Wild-type fitness is 1.'],
 ['Fitness','LD starts day −5. CAR-associated fitness rises before Tmax and decays after Tmax. Shared kinetics are fixed in each scenario.'],
 ['Double hit','Second transforming alteration, aggregated across double-hit cells. Double-hit fitness is fM + sD. Baseline double hits may exist.'],
 ['Fitting','One s₂ per patient/gene from first and last VAF only; pseudo-Huber log-odds loss. All eligible variants share the patient/gene estimate.'],
 ['VAF mapping','Founder-marked cell fraction = 2 × VAF fraction under a diploid heterozygous assumption. Zero, missing, or ≥50% VAF is not fitted.'],
 ['Variant dependence','Do not add variant VAFs or combine gene-specific risks without clone phasing. Variants can belong to the same clone.'],
 ['Uncertainty','60,000 patient-cluster bootstrap draws, retaining all gene memberships together. Fixed kinetics and measurement uncertainty are not included.'],
 ['Multiplicity','Bonferroni-adjusted percentile intervals cover the 15 planned contrasts approximately; the small groups limit inference. No sequencing likelihood is claimed.'],
 ['Risk method','Hybrid deterministic single-hit background plus an adaptive absorbing double-hit birth/death forward solver. This approximates full three-type stochastic Moran dynamics.'],
 ['Numerical precision','Solver rtol = 1e−7 and per-state atol = 1e−12. Probabilities below atol lack a relative-accuracy guarantee. The zero at μ=0 and D0=0 is exact.'],
 ['Risk conditioning','Patient risk begins at the first endpoint, conditional on no earlier threshold crossing; s₂ uses the later VAF retrospectively and remains fixed across scenarios.'],
 ['Endpoint diagnostics','Coupled mean-field endpoint checks in the patient risk table reveal when an assumed transforming scenario changes the fitted founder burden. No refit is implied.'],
 ['Clinical comparison','Source status is descriptive. Event-free follow-up is unavailable; the clinical tMN label is not calibrated to the modeled 1% threshold.'],
 ['Refresh','Run the accompanying extraction, fitting and risk scripts, then rebuild this static workbook from their result tables.'],
 ['Assumption provenance','Reference kinetics and sensitivity/risk grids were approved by the user on 2026-11-05. They are modeling assumptions, not measured progenitor parameters.']
];
settings.getRangeByIndexes(40,0,notes.length+1,2).values=[['Definition','Details'],...notes];
settings.getRange(`A41:B${notes.length+41}`).format.font={name:'Arial',size:10,color:ink};
settings.getRange('A41:D41').format={fill:pale,font:{bold:true,color:navy}};
settings.getRange(`B42:B${notes.length+41}`).format.wrapText=false;
payload.tables.find(t=>t.name==='Model Settings').rows.forEach((row,i)=>{settings.getRange(`B${i+2}`).setNumberFormat(row[0].startsWith('mu_')?'0.00E+00':Number.isInteger(row[1])?'#,##0':'0.0000');});
wb.recalculate();
await fs.writeFile(path.join(base,'workbook_manifest.json'),JSON.stringify(manifest,null,2));
const inspect=await wb.inspect({kind:'table',range:'Summary!A6:H12',include:'values,formulas',tableMaxRows:8,tableMaxCols:8,maxChars:5000});
await fs.writeFile(path.join(base,'summary_inspection.ndjson'),inspect.ndjson);
const errors=await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:20},maxChars:2000});
await fs.writeFile(path.join(base,'formula_error_scan.ndjson'),errors.ndjson);
console.log('Inspected formulas');
await fs.mkdir(path.join(base,'previews'),{recursive:true});
for(const t of payload.tables){
 const range=t.name==='Summary'?'A1:H28':t.name==='Model Settings'?'A1:D8':'A1:H8';
 const png=await wb.render({sheetName:t.name,range,scale:1.2,format:'png'});
 await fs.writeFile(path.join(base,'previews',t.name.replaceAll(' ','_')+'.png'),new Uint8Array(await png.arrayBuffer()));
 console.log(`Rendered ${t.name}`);
}
if(!previewOnly&&!process.argv.includes('--defer-export')){
 await fs.mkdir(path.dirname(output),{recursive:true});
 const xlsx=await SpreadsheetFile.exportXlsx(wb);await xlsx.save(output);
 console.log(`EXPORTED ${output}`);
}else console.log('Preview check completed; no workbook exported.');
