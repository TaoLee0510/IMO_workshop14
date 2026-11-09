"""Independent read-only reconciliation of prepared endpoints against approved sheets."""
import csv
import json
import math
from collections import defaultdict
from pathlib import Path
import openpyxl

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT.parent

def main():
    rows = list(csv.DictReader((ROOT / 'data/all_endpoints.csv').open()))
    w9 = openpyxl.load_workbook(SOURCE / '41588_2026_2795_MOESM9_ESM.xlsx', read_only=True, data_only=True)
    raw = list(w9['Gene_Level_Expansion_Data'].values)
    w9.close()
    w7 = openpyxl.load_workbook(SOURCE / '41588_2026_2795_MOESM7_ESM.xlsx', read_only=True, data_only=True)
    status = list(w7['tMN_VAF'].values)
    clinical = list(w7['Auto_CAR_tMN'].values)
    w7.close()
    positives = defaultdict(list)
    for i, r in enumerate(status[1:], 2):
        if r[-1] == 'tMN': positives[r[0]].append(i)
    checked = 0
    for r in rows:
        source = raw[int(r['source_excel_row']) - 1]
        assert r['patient_id'] == source[0] and r['gene'] == source[1] == source[8]
        assert r['variant_id'] == f'{source[9]}:{source[10]}:{source[11]}'
        for key, pos in [('first_day', 2), ('last_day', 3), ('first_vaf_percent', 4), ('last_vaf_percent', 5)]:
            value = source[pos]
            if isinstance(value, (int, float)):
                assert float(r[key]) == value, (key, r['source_excel_row'])
            else: assert r[key] == '', (key, r['source_excel_row'])
        for key, pos in [('first_vaf_source_value', 4), ('last_vaf_source_value', 5)]:
            value = source[pos]
            assert (float(r[key]) == value if isinstance(value, (int, float)) else r[key] == str(value))
        assert int(r['tmn_source_flag']) == int(source[0] in positives)
        assert r['tmn_source_excel_rows'] == ';'.join(map(str, positives.get(source[0], [])))
        eligible = all(isinstance(source[k], (int, float)) for k in (2,3,4,5)) and source[2] < source[3] and 0 < source[4] < 50 and 0 < source[5] < 50
        assert (r['eligible_endpoint_pair'] == 'True') == eligible
        checked += 1
    assert len(raw)-1 == checked == 171
    inventory = list(csv.DictReader((ROOT / 'data/gene_eligibility.csv').open()))
    for g in inventory:
        rr = [r for r in rows if r['gene'] == g['gene']]
        ee = [r for r in rr if r['eligible_endpoint_pair'] == 'True']
        assert int(g['all_variant_records']) == len(rr)
        assert int(g['all_patient_gene_groups']) == len({r['patient_id'] for r in rr})
        assert int(g['eligible_variant_records']) == len(ee)
        assert int(g['eligible_patient_gene_groups']) == len({r['patient_id'] for r in ee})
        assert (g['retained_by_gt2_filter'] == 'True') == (len({r['patient_id'] for r in ee}) > 2)
    outcomes = json.loads((ROOT / 'data/outcome_evidence.json').read_text())
    for o in outcomes:
        expected = [{'source_excel_row': i, 'cohort': r[1], 'time_cell_therapy_to_tMN': r[2]} for i,r in enumerate(clinical[1:],2) if r[0] == o['patient_id']]
        assert o['clinical_records'] == expected
        assert o['tmn_source_flag'] == int(o['patient_id'] in positives)
    result = {'passed': True, 'source_endpoint_rows_reconciled': checked, 'gene_inventory_rows_reconciled':len(inventory),
              'outcome_patients_reconciled':len(outcomes), 'source_sheets_read':['MOESM9:Gene_Level_Expansion_Data','MOESM7:tMN_VAF','MOESM7:Auto_CAR_tMN'],
              'longitudinal_observations_read':False, 'source_NA_values_preserved':sum(r['first_vaf_source_value']=='NA' for r in rows)}
    (ROOT/'workbook_support/source_reconciliation.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result))

if __name__ == '__main__': main()
