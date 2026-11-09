"""Read-only full-cell XLSX comparison against typed CSV/JSON export payload."""
import json
import math
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path
import openpyxl

ROOT=Path(__file__).resolve().parents[1]
PATH=ROOT.parent/'outputs/01a1ac40-905d-7660-b1ab-5500bfe0661f/tmn_all_genes_gt2_results.xlsx'

def main():
    payload=json.loads((ROOT/'workbook_support/workbook_payload.json').read_text())
    workbook=openpyxl.load_workbook(PATH,read_only=True,data_only=True)
    differences=[]; counts={};numeric=0;text=0;max_abs=0;max_rel=0
    for t in payload['tables']:
        sheet=workbook[t['name']]; start=7 if t['name']=='Summary' else 2
        rows=sheet.iter_rows(min_row=start,max_row=start+len(t['rows'])-1,max_col=len(t['columns']),values_only=True)
        for i,(expected,actual) in enumerate(zip(t['rows'],rows),start):
            for j,(a,b) in enumerate(zip(expected,actual),1):
                if isinstance(a,(int,float)) and not isinstance(a,bool):
                    numeric+=1
                    if not isinstance(b,(int,float)) or isinstance(b,bool): differences.append([t['name'],i,j,a,b,'numeric_type'])
                    else:
                        absolute=abs(a-b);relative=absolute/max(abs(a),1e-300)
                        max_abs=max(max_abs,absolute);max_rel=max(max_rel,relative)
                        if not math.isclose(a,b,rel_tol=2e-15,abs_tol=1e-300):differences.append([t['name'],i,j,a,b,'numeric_precision'])
                else:
                    text+=1
                    if a!=b or (isinstance(a,bool) and not isinstance(b,bool)):differences.append([t['name'],i,j,a,b,'text_blank_boolean'])
        counts[t['name']]=len(t['rows'])
    summary=[workbook['Summary'][f'D{i}'].value for i in range(15,21)]
    assert summary==[171,108,83,110,28,31],summary
    assert workbook['Endpoints']['F2'].value==0
    na=[row for row in workbook['Endpoints'].iter_rows(min_row=2,max_row=172,values_only=True) if row[0]=='Patient_103' and row[1]=='ASXL1'][0]
    assert na[5] is None and na[7]=='NA'
    workbook.close()
    ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    panes=[]
    with zipfile.ZipFile(PATH) as z:
        for i,t in enumerate(payload['tables'],1):
            xml=ET.fromstring(z.read(f'xl/worksheets/sheet{i}.xml'))
            pane=xml.find('s:sheetViews/s:sheetView/s:pane',ns)
            panes.append({'sheet':t['name'],'pane':pane.attrib if pane is not None else None})
            if t['name']!='Summary':
                assert pane is not None and pane.get('state')=='frozen'
                assert float(pane.get('ySplit','0'))>=1,(t['name'],pane.attrib)
                assert float(pane.get('xSplit','0'))>=1,(t['name'],pane.attrib)
        tables=[n for n in z.namelist() if n.startswith('xl/tables/table') and n.endswith('.xml')]
        assert len(tables)==len(payload['tables'])
        for file in tables:assert ET.fromstring(z.read(file)).find('s:autoFilter',ns) is not None
    result={'passed':not differences,'sheet_row_counts':counts,'numeric_cells_checked':numeric,'other_cells_checked':text,'max_absolute_numeric_difference':max_abs,'max_relative_numeric_difference':max_rel,'differences':differences[:20],'difference_count':len(differences),'summary_formulas_verified':summary,'tables_with_autofilters':len(tables),'frozen_panes':panes,'NA_and_zero_preserved':True}
    (ROOT/'workbook_support/workbook_verification.json').write_text(json.dumps(result,indent=2))
    print(json.dumps({k:v for k,v in result.items() if k not in ('frozen_panes','sheet_row_counts')}))
    assert not differences,'Export changed typed source values; see verification report.'

if __name__=='__main__':main()
