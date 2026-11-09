import fs from 'node:fs/promises';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import {Workbook,SpreadsheetFile} from '@oai/artifact-tool';
const root=path.dirname(fileURLToPath(import.meta.url));
const data=JSON.parse(await fs.readFile(path.join(root,'results/workbook_payload.json'),'utf8'));
const output=path.join(root,'outputs/20261105_patient_analysis');await fs.mkdir(output,{recursive:true});
const wb=Workbook.create();const sheets=[];
function table(name,rows,keys,start=0){
 const sh=wb.worksheets.add(name);sheets.push(name);sh.showGridLines=false;
 const matrix=[keys,...rows.map(r=>keys.map(k=>r[k]??null))];
 const rg=sh.getRangeByIndexes(start,0,matrix.length,keys.length);rg.values=matrix;
 rg.format.font={name:'Arial',size:10};rg.format.rowHeight=23;rg.format.columnWidth=19;
 const header=sh.getRangeByIndexes(start,0,1,keys.length);
 header.format={fill:'#193f5e',font:{name:'Arial',size:10,bold:true,color:'#FFFFFF'},wrapText:true,rowHeight:48};
 sh.getRangeByIndexes(start,0,matrix.length,1).format.columnWidth=23;
 sh.freezePanes.freezeRows(start+1);sh.freezePanes.freezeColumns(1);
 keys.forEach((k,j)=>{
  const col=sh.getRangeByIndexes(start+1,j,rows.length,1);
  if(rows.some(r=>typeof r[k]==='number'))col.setNumberFormat('0.0000');
  if(/probability/.test(k)&&!k.startsWith('log10'))col.setNumberFormat('0.000E+00');
  if(/(^|_)VAF$|fraction/.test(k)&&!/count/.test(k))col.setNumberFormat('0.0000%');
  if(/initial_[WCD]|_cells|replicates|_count/.test(k))col.setNumberFormat('#,##0.00');
  if(/initial_[WCD]|replicates|_count/.test(k))col.setNumberFormat('#,##0');
  if(k==='day'||k==='baseline_day'||k==='outcome_day')col.setNumberFormat('0');
  if(k==='source'||k==='definition')col.format.columnWidth=65;
  if(k==='variant_id')col.format.columnWidth=32;
 });
 return sh;
}
const riskKeys=['patient_id','gene','baseline_VAF','year1_probability','year2_probability','year3_probability','year4_probability','outcome'];
const risk=table('Patient risks',data.risks,riskKeys,3);
risk.getRange('A2').values=[['Patient-specific tMN model probabilities']];risk.getRange('A2').format.font={name:'Arial',size:14,bold:true};
risk.getRange('A21').values=[['Forecast origin']];risk.getRange('B21').values=[['CAR infusion day 0']];
risk.getRange('A22').values=[['Conditioning']];risk.getRange('B22').values=[['Patient baseline and event-free day 0']];
risk.getRange('B21:H22').format.columnWidth=24;risk.getRange('B22').format.columnWidth=47;
table('VAF observations',data.observations,['patient_id','gene','day','observed_VAF','mean_VAF','lower_VAF','upper_VAF','PIT','predictive_residual','approximate_z_residual']);
const ll=table('Likelihoods',data.likelihoods,['patient_id','gene','outcome','VAF_nll','clinical_given_VAF_nll','joint_nll']);
ll.getRange('A16').values=[['Total']];ll.getRange('D16:F16').formulas=[['=SUM(D2:D15)','=SUM(E2:E15)','=SUM(F2:F15)']];ll.getRange('A16:F16').format.font={bold:true};
table('Risk curves',data.risk_curves,['patient_id','day','probability','log10_probability','numeric_underflow']);
table('Simulation summary',data.simulations,['patient_id','replicates','initial_C','initial_D','mean_C_year4','mean_D_year4','fraction_ever_D','MC_probability_year1','MC_probability_year4']);
table('Yearly trajectories',data.yearly_trajectories,['patient_id','day','mean_W_cells','mean_C_cells','mean_D_cells','lower_C_cells','upper_C_cells','lower_D_cells','upper_D_cells']);
table('Initial states',data.initial_states,['patient_id','gene','sample_type','baseline_day','baseline_VAF','initial_W','initial_C','initial_D','rounding_C_cells','outcome','outcome_day','followup_VAF_count','variant_id','source']);
const parameters=table('Parameters',data.parameters,['parameter','value','definition']);
parameters.getRange('B2:B11').setNumberFormat('0.000000');
data.parameters.forEach((r,i)=>parameters.getRange('B'+(i+2)).setNumberFormat(r.parameter==='mu'?'0.000E+00':['N','D_init'].includes(r.parameter)?'#,##0':'0.000000'));
parameters.getRange('A14:D16').values=[
 ['Input source','MOESM9','Filtered mutation data; baseline-selected founder marker and longitudinal VAF','Source rows retained in cohort input'],
 ['Clinical source','MOESM7','Auto_CAR_tMN rows 33, 37, 31; events at 127 and 589 days','Patient 25: listed first event day 1561; censor at day 1460'],
 ['Risk definition','D > 1%','Strict aggregate threshold D >= 1001; conditional model probabilities','Outcomes unknown for 11 patients']
];
parameters.getRange('A14:D16').format.rowHeight=48;parameters.getRange('B14:D16').format.wrapText=true;parameters.getRange('C1:D20').format.columnWidth=64;
data.notes.forEach((note,i)=>{parameters.getRange('A'+(19+i)).values=[['Model note']];parameters.getRange('C'+(19+i)).values=[[note]];parameters.getRange('C'+(19+i)).format.wrapText=true;parameters.getRange('A'+(19+i)+':D'+(19+i)).format.rowHeight=43;});
wb.recalculate();
await fs.writeFile(path.join(root,'results/workbook_inspection.json'),JSON.stringify(await wb.inspect({kind:'table',range:'Patient risks!A4:H18',include:'values,formulas',tableMaxRows:15,tableMaxCols:8})));
await fs.writeFile(path.join(root,'results/workbook_formula_scan.json'),JSON.stringify(await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#NUM!',options:{useRegex:true,maxResults:30}})));
for(const name of sheets){
 if(process.env.PREVIEW_SHEETS&&!process.env.PREVIEW_SHEETS.split(',').includes(name))continue;
 const range=name==='Patient risks'?'A2:H18':name==='Parameters'?'A1:D22':'A1:H12';
 const preview=await wb.render({sheetName:name,range,scale:1.2,format:'png'});
 await fs.writeFile(path.join(output,'preview-'+name.replaceAll(' ','-')+'.png'),new Uint8Array(await preview.arrayBuffer()));
}
const xlsx=await SpreadsheetFile.exportXlsx(wb);await xlsx.save(path.join(output,'Patient_specific_Moran_results.xlsx'));console.log('Workbook exported.');
