import fs from 'node:fs/promises';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import {Workbook,SpreadsheetFile} from '@oai/artifact-tool';
const root=path.dirname(fileURLToPath(import.meta.url));
const data=JSON.parse(await fs.readFile(path.join(root,'results/workbook_payload.json'),'utf8'));
const output=path.join(root,'outputs/20261105_mu001_sensitivity');await fs.mkdir(output,{recursive:true});
const wb=Workbook.create();const summary=wb.worksheets.add('Summary');summary.showGridLines=false;const names=['Summary'];
function table(name,rows,keys){
 const sh=wb.worksheets.add(name);names.push(name);sh.showGridLines=false;
 const values=[keys,...rows.map(r=>keys.map(k=>r[k]??null))];const rg=sh.getRangeByIndexes(0,0,values.length,keys.length);
 rg.values=values;rg.format.font={name:'Arial',size:10};rg.format.rowHeight=22;rg.format.columnWidth=19;
 sh.getRangeByIndexes(0,0,1,keys.length).format={fill:'#193f5e',font:{name:'Arial',size:10,bold:true,color:'#FFFFFF'},wrapText:true,rowHeight:45};
 sh.getRangeByIndexes(0,0,values.length,1).format.columnWidth=27;sh.freezePanes.freezeRows(1);sh.freezePanes.freezeColumns(1);
 keys.forEach((k,j)=>{
  const col=sh.getRangeByIndexes(1,j,rows.length,1);
  if(rows.some(r=>typeof r[k]==='number'))col.setNumberFormat('0.0000');
  if(/probability/.test(k))col.setNumberFormat('0.000E+00');
  if(/(^|_)VAF$|fraction/.test(k))col.setNumberFormat('0.0000%');
  if(k==='mu'||/refinement|difference/.test(k))col.setNumberFormat('0.000E+00');
  if(/D_init|initial_[WCD]_integer|^initial_[WCD]$|replicates|evaluation_count|error_count/.test(k))col.setNumberFormat('#,##0');
  if(k==='day'||k==='baseline_day')col.setNumberFormat('0');
  if(k==='definition')col.format.columnWidth=66;
 });return sh;
}
table('Joint fits',data.fits.map(r=>({...r,boundary_optimum:r.boundary_optimum==null?null:r.boundary_optimum?'At bound':'Interior'})),['scenario','status','mu','D_init','alpha','theta_days','joint_nll','VAF_nll','clinical_given_VAF_nll','b','s_D','peak_day','boundary_optimum','evaluation_count','cap_refinement_loglik_difference','tolerance_refinement_loglik_difference','full_N_loglik_difference']);
table('Patient risks',data.risks,['scenario','patient_id','gene','baseline_VAF','baseline_day','initial_D','mu','year1_probability','year2_probability','year3_probability','year4_probability']);
table('Reference VAF',data.observations,['patient_id','gene','day','observed_VAF','mean_VAF','lower_VAF','upper_VAF','PIT','predictive_residual','approximate_z_residual']);
table('Reference sims',data.simulations,['patient_id','replicates','initial_C','initial_D','mean_C_year4','mean_D_year4','fraction_ever_D','MC_probability_year1','MC_probability_year4']);
table('Reference yearly',data.yearly,['patient_id','day','mean_W_cells','mean_C_cells','mean_D_cells','lower_C_cells','upper_C_cells','lower_D_cells','upper_D_cells','hybrid_mean_C_cells','hybrid_mean_D_cells','hybrid_tMN_probability','MC_tMN_probability','MC_probability_lower_95','MC_probability_upper_95']);
table('Initial states',data.initial_states,['patient_id','gene','baseline_day','baseline_VAF','D_init','initial_C_continuous','initial_C_integer','initial_W_integer','D_fraction_of_founder','sample_type']);
const parameters=table('Parameters',data.parameters,['parameter','value','definition']);
parameters.getRange('B2:B11').setNumberFormat('0.000000');
data.parameters.forEach((r,i)=>parameters.getRange('B'+(i+2)).setNumberFormat(r.parameter==='mu'?'0.000E+00':['N','D_init','K_cells','diagnosis_delay_days'].includes(r.parameter)?'#,##0':'0.000000'));
parameters.getRange('A17:D20').values=[
 ['VAF source','MOESM9','Filtered mutation data; source observation rows retained in the frozen cohort','One baseline-selected founder marker per patient'],
 ['Clinical source','MOESM7','Auto_CAR_tMN rows 33, 37, 31: events at 127 and 589 days','Patient 25 censored at 1460 from its listed first event at 1561'],
 ['Turnover precedent','Heyde et al. (2021)','https://pmc.ncbi.nlm.nih.gov/articles/PMC8109274/','Slow range 1/280 to 1/28 per day; 1/7 is a fast-turnover stress case'],
 ['Scenario interpretation','Conditional forecasts','Same cohort and VAF/tMN likelihood; mu and initial D are fixed in each scenario','Unknown clinical outcomes are not negatives']
 ];
parameters.getRange('C1:D20').format.columnWidth=66;parameters.getRange('A17:D20').format.rowHeight=48;parameters.getRange('B17:D20').format.wrapText=true;
summary.getRange('A2').values=[['Mu=0.01 joint refit and sensitivity']];summary.getRange('A2').format.font={name:'Arial',size:14,bold:true};
summary.getRange('A4:C4').values=[['Metric','Value','Meaning']];summary.getRange('A4:C4').format={fill:'#193f5e',font:{name:'Arial',size:10,bold:true,color:'#FFFFFF'}};
const metrics=[['Fitted alpha',null,'Shared Gamma shape'],['Fitted theta (days)',null,'Shared Gamma scale'],['Reference mu',null,'C-to-D probability per C division'],['Reference initial D',null,'Cells at measured baseline'],['Joint NLL',null,'VAF plus clinical-given-VAF NLL'],['VAF NLL',null,'Marginal VAF contribution'],['Clinical given VAF NLL',null,'Clinical factor on shared latent states'],['Total scenarios',79,'72 mu-D grid plus 7 one-at-a-time cases'],['Compatible fitted scenarios',null,'Alpha and theta re-estimated in every compatible case'],['Structurally incompatible',null,'Zero mutation and zero initial D'],['Follow-up VAF observations',null,'Same 14-patient cohort']];
summary.getRange('A5:C15').values=metrics;summary.getRange('A4:C15').format.font={name:'Arial',size:10};summary.getRange('A4:C15').format.rowHeight=25;
summary.getRange('A4:A34').format.columnWidth=31;summary.getRange('B4:B34').format.columnWidth=20;summary.getRange('C4:C34').format.columnWidth=54;
summary.getRange('B5:B11').formulas=[["='Joint fits'!E2"],["='Joint fits'!F2"],["='Joint fits'!C2"],["='Joint fits'!D2"],["='Joint fits'!G2"],["='Joint fits'!H2"],["='Joint fits'!I2"]];
summary.getRange('B5:B11').setNumberFormat('0.000000');summary.getRange('B8').setNumberFormat('0');
summary.getRange('B13').formulas=[['=COUNTIFS(\'Joint fits\'!B2:B80,"fitted")']];
summary.getRange('B14').formulas=[['=COUNTIFS(\'Joint fits\'!B2:B80,"structurally_incompatible")']];
summary.getRange('B15').formulas=[["=COUNT('Reference VAF'!C2:C43)"]];
summary.getRange('A17:C18').values=[['95% VAF interval coverage',data.diagnostics.predictive_interval_coverage_fraction,'Observed coverage: 22 of 42 follow-up VAF measurements'],['Fits at a search boundary',null,'Conditional estimates; bounds do not identify a physiological timescale']];
summary.getRange('B17').setNumberFormat('0.0%');summary.getRange('B18').formulas=[['=COUNTIF(\'Joint fits\'!M2:M80,"At bound")']];
summary.getRange('C17:C18').format.wrapText=true;summary.getRange('A17:C18').format.rowHeight=30;
summary.getRange('A20:D20').values=[['Patient','Gene','Reference one-year probability','Reference four-year probability']];
summary.getRange('A20:D20').format={fill:'#193f5e',font:{name:'Arial',size:10,bold:true,color:'#FFFFFF'},wrapText:true,rowHeight:34};
summary.getRange('D20:D34').format.columnWidth=29;summary.getRange('A21:D34').format.rowHeight=23;
const links=data.reference_risks.map((r,i)=>["='Patient risks'!B"+(i+2),"='Patient risks'!C"+(i+2),"='Patient risks'!H"+(i+2),"='Patient risks'!K"+(i+2)]);
summary.getRange('A21:D34').formulas=links;summary.getRange('C21:D34').setNumberFormat('0.000E+00');
wb.recalculate();
await fs.writeFile(path.join(root,'results/workbook_inspection.json'),JSON.stringify(await wb.inspect({kind:'table',range:'Summary!A4:C15',include:'values,formulas',tableMaxRows:12,tableMaxCols:3})));
await fs.writeFile(path.join(root,'results/workbook_formula_scan.json'),JSON.stringify(await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:30}})));
for(const name of names){
 if(process.env.PREVIEW_SHEETS&&!process.env.PREVIEW_SHEETS.split(',').includes(name))continue;
 const range=({'Summary':'A2:D34','Parameters':'A1:D20','Patient risks':'A1:K12','Joint fits':'A1:Q12','Reference yearly':'A1:O12','Initial states':'A1:J12','Reference VAF':'A1:J12'})[name]||'A1:I12';
 const preview=await wb.render({sheetName:name,range,scale:1.2,format:'png'});
 await fs.writeFile(path.join(output,'preview-'+name.replaceAll(' ','-')+'.png'),new Uint8Array(await preview.arrayBuffer()));
}
const outputFile=await SpreadsheetFile.exportXlsx(wb);await outputFile.save(path.join(output,'Mu001_joint_sensitivity_results.xlsx'));console.log('Workbook exported.');
