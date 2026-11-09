import fs from 'node:fs/promises';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import {Workbook,SpreadsheetFile} from '@oai/artifact-tool';
const root=path.dirname(fileURLToPath(import.meta.url));
const data=JSON.parse(await fs.readFile(path.join(root,'results/workbook_payload.json'),'utf8'));
const output=path.join(root,'outputs/20261105_fixed_area');
await fs.mkdir(output,{recursive:true});
const wb=Workbook.create();
function table(name,rows,columns){
 const sh=wb.worksheets.add(name);sh.showGridLines=false;
 const keys=columns??Object.keys(rows[0]);
 const values=[keys,...rows.map(r=>keys.map(k=>r[k]??null))];
 const rg=sh.getRangeByIndexes(0,0,values.length,keys.length);rg.values=values;
 rg.format.font={name:'Arial',size:10};rg.format.rowHeight=21;rg.format.columnWidth=17;
 const head=sh.getRangeByIndexes(0,0,1,keys.length);
 head.format={fill:'#193f5e',font:{name:'Arial',size:10,bold:true,color:'#FFFFFF'},wrapText:true,rowHeight:42};
 sh.freezePanes.freezeRows(1);sh.freezePanes.freezeColumns(1);
 for(let c=0;c<keys.length;c++){
  if(rows.some(r=>typeof r[keys[c]]==='number'))sh.getRangeByIndexes(1,c,rows.length,1).setNumberFormat('0.0000');
  if(keys[c].includes('probability'))sh.getRangeByIndexes(1,c,rows.length,1).setNumberFormat('0.000E+00');
  if(keys[c].includes('fraction')||keys[c].includes('MC_risk'))sh.getRangeByIndexes(1,c,rows.length,1).setNumberFormat('0.0000%');
  if(keys[c]==='mu'||keys[c].includes('refinement')||keys[c].includes('roundoff')||keys[c].includes('_SE'))sh.getRangeByIndexes(1,c,rows.length,1).setNumberFormat('0.000E+00');
  if(['N','D_init','replicates','patients','VAF_observations','followup_VAF_count','cap','nfev'].includes(keys[c]))sh.getRangeByIndexes(1,c,rows.length,1).setNumberFormat('#,##0');
  if(keys[c]==='b')sh.getRangeByIndexes(1,c,rows.length,1).setNumberFormat('0.000000');
  if(keys[c]==='baseline_vaf')sh.getRangeByIndexes(1,c,rows.length,1).setNumberFormat('0.0000%');
  if(keys[c]==='mean_D_four_year_fraction')sh.getRangeByIndexes(1,c,rows.length,1).setNumberFormat('0.000000%');
 }
 sh.getRangeByIndexes(1,0,rows.length,1).format.columnWidth=30;
 return sh;
}
const summary=wb.worksheets.add('Results');summary.showGridLines=false;
summary.getRange('A2').values=[['Fixed-area gamma Moran study']];summary.getRange('A2').format.font={name:'Arial',size:15,bold:true};
summary.getRange('A4:C4').values=[['Metric','Value','Meaning']];summary.getRange('A4:C4').format={fill:'#193f5e',font:{name:'Arial',bold:true,color:'#FFFFFF'}};
const ref=data.reference;
const summaryRows=[['Clinical events used',2,'Exact event times: 127 and 589 days'],['Four-year censor records',1,'First recorded tMN at day 1561'],['Unknown clinical outcomes',11,'Not coded as negatives'],['Patients',14,'One founder marker per patient'],['Follow-up VAF observations',42,'Baseline conditions are separate'],['Gamma advantage area (days)',1,'Fixed area; peak is free to change with shape/scale'],['Fitted alpha',ref.alpha,'Shared gamma shape'],['Fitted theta (days)',ref.theta_days,'Shared gamma scale'],['Peak day',ref.peak_day,'Derived, not independently fitted'],['Joint NLL',ref.joint_nll,'VAF plus clinical-given-VAF NLL'],['VAF NLL',ref.VAF_nll,'Marginal VAF contribution'],['Clinical-given-VAF NLL',ref.clinical_given_VAF_nll,'Clinical contribution on the shared latent process']];
summary.getRange('A5:C16').values=summaryRows;
summary.getRange('A4:C16').format.font={name:'Arial',size:10};summary.getRange('A4:C16').format.rowHeight=24;
summary.getRange('A4:A16').format.columnWidth=34;summary.getRange('B4:B16').format.columnWidth=19;summary.getRange('C4:C16').format.columnWidth=67;
summary.getRange('B5:B16').setNumberFormat('0.0000');
summary.getRange('B5:B9').setNumberFormat('#,##0');summary.getRange('B10').setNumberFormat('0.0');
summary.getRange('A19').values=[['Interpretation']];summary.getRange('A19').format.font={bold:true};
data.notes.forEach((note,i)=>{summary.getRange(`A${20+i}`).values=[['Scope / interpretation']];summary.getRange(`C${20+i}`).values=[[note]];summary.getRange(`C${20+i}`).format.wrapText=true;summary.getRange(`A${20+i}:C${20+i}`).format.rowHeight=38;});
const fit=table('Joint fits',data.fit_rows);
table('Likelihoods',data.likelihood_rows);
table('Four-year risks',data.risk_rows);
table('Simulations',data.simulation_rows);
table('Patients',data.patient_rows);
const parameters=table('Parameters',data.parameter_rows);parameters.getRange('D1:D12').format.columnWidth=65;parameters.getRange('E1:E12').format.columnWidth=80;
parameters.getRange('D2:E10').format.wrapText=true;parameters.getRange('A2:E10').format.rowHeight=36;
data.parameter_rows.forEach((r,i)=>parameters.getRange(`B${i+2}`).setNumberFormat(r.parameter==='mu'?'0.00E+00':['N','D_init'].includes(r.parameter)?'#,##0':'0.000000'));
parameters.getRange('A14:E15').values=[['Clinical source','41588_2026_2795_MOESM7_ESM.xlsx',null,'Auto_CAR_tMN: rows 33, 37 and 31. Patient 25 is censored at day 1460 using its listed first event at day 1561.',null],['VAF source','41588_2026_2795_MOESM9_ESM.xlsx',null,'Filtered_Mutation_Data; index-marker extraction and exact longitudinal row references are retained in the study data files.',null]];
parameters.getRange('A14:E15').format.rowHeight=48;parameters.getRange('B14:D15').format.wrapText=true;
if(data.profile_rows.length)table('Profiles',data.profile_rows);
if(data.crosscheck_rows.length)table('MC comparison',data.crosscheck_rows);
for(const [name,rows] of [['Joint fits',data.fit_rows],['Four-year risks',data.risk_rows],['Simulations',data.simulation_rows],['Patients',data.patient_rows]]){
 const sh=wb.worksheets.getItem(name);
 const keys=Object.keys(rows[0]);
 for(let c=0;c<keys.length;c++)if(['scenario','status','name','label','patient_id','variant_id','source_evidence','clinical_source_sheet'].includes(keys[c]))sh.getRangeByIndexes(0,c,rows.length+1,1).format.columnWidth=30;
}
wb.recalculate();
await fs.writeFile(path.join(root,'results/workbook_inspection.json'),JSON.stringify(await wb.inspect({kind:'table',range:'Results!A4:C16',include:'values,formulas',tableMaxRows:13,tableMaxCols:3})));
const errors=await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#NUM!',options:{useRegex:true,maxResults:30}});
await fs.writeFile(path.join(root,'results/workbook_formula_scan.json'),JSON.stringify(errors));
for(const name of ['Results','Joint fits','Likelihoods','Four-year risks','Simulations','Patients','Parameters',...(data.profile_rows.length?['Profiles']:[]),...(data.crosscheck_rows.length?['MC comparison']:[])]){
 const preview=await wb.render({sheetName:name,range:name==='Results'?'A2:C26':name==='Parameters'?'A1:E15':'A1:H12',scale:1.3,format:'png'});
 await fs.writeFile(path.join(output,`preview-${name.replaceAll(' ','-')}.png`),new Uint8Array(await preview.arrayBuffer()));
}
const xlsx=await SpreadsheetFile.exportXlsx(wb);
await xlsx.save(path.join(output,'Fixed_area_joint_Moran_study.xlsx'));
console.log('Workbook exported.');
