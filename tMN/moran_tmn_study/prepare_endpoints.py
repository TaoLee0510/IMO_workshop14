"""Extract only the approved first/last endpoints and source outcome evidence."""
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
import openpyxl

ROOT = Path(__file__).resolve().parent
DATA = ROOT / 'data'
DATA.mkdir(exist_ok=True)
SOURCE = ROOT.parent
GENES = ('DNMT3A', 'TP53', 'PPM1D')

def main():
    f7 = SOURCE / '41588_2026_2795_MOESM7_ESM.xlsx'
    f9 = SOURCE / '41588_2026_2795_MOESM9_ESM.xlsx'
    w7 = openpyxl.load_workbook(f7, read_only=True, data_only=True)
    vaf_rows = list(w7['tMN_VAF'].values)
    assert vaf_rows[0][0] == 'Paper_ID' and vaf_rows[0][-1] == 'Timpoint'
    positive = {r[0] for r in vaf_rows[1:] if r[-1] == 'tMN'}
    clinical = defaultdict(list)
    for excel_row, row in enumerate(w7['Auto_CAR_tMN'].values, 1):
        if excel_row > 1:
            clinical[row[0]].append({'source_excel_row': excel_row, 'cohort': row[1], 'time_cell_therapy_to_tMN': row[2]})
    w7.close()
    w9 = openpyxl.load_workbook(f9, read_only=True, data_only=True)
    # Do not read Filtered_Mutation_Data or any intermediate VAF measurements.
    records = []
    for excel_row, row in enumerate(w9['Gene_Level_Expansion_Data'].values, 1):
        if excel_row == 1 or row[1] not in GENES:
            continue
        assert row[1] == row[8]
        assert row[2] < row[3]
        flags = []
        if row[4] == 0:
            flags.append('baseline_zero_requires_detection_model')
        if row[5] == 0:
            flags.append('last_zero_requires_detection_model')
        if max(row[4:6]) >= 50:
            flags.append('diploid_heterozygous_mapping_not_valid_for_fit')
        if row[5] < row[4]:
            flags.append('decreasing_endpoint_nonnegative_selection_mismatch')
        if row[2] > -5:
            flags.append('baseline_after_LD_start')
        if clinical.get(row[0]) and row[0] not in positive:
            flags.append('clinical_tMN_record_without_source_flag')
        clinical_days = [r['time_cell_therapy_to_tMN'] for r in clinical[row[0]] if r['cohort'] == 'CAR']
        if len(clinical_days) > 1:
            flags.append('multiple_clinical_time_origins_require_reconciliation')
        if any(isinstance(d, (int, float)) and row[3] >= d for d in clinical_days):
            flags.append('last_VAF_at_or_after_listed_diagnosis_review_required')
        if row[0] in positive and not clinical[row[0]]:
            flags.append('source_positive_without_clinical_timing_record')
        records.append({
            'patient_id': row[0], 'gene': row[1],
            'variant_id': f'{row[9]}:{row[10]}:{row[11]}',
            'first_day': row[2], 'last_day': row[3],
            'first_vaf_percent': row[4], 'last_vaf_percent': row[5],
            'source_excel_row': excel_row,
            'tmn_source_flag': int(row[0] in positive),
            'quality_flags': ';'.join(flags),
        })
    w9.close()
    assert len(records) == 130
    assert len({(r['patient_id'], r['gene'], r['variant_id']) for r in records}) == 130
    with (DATA / 'endpoints.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    ids = {r['patient_id'] for r in records}
    outcomes = []
    for pid in sorted(ids, key=lambda p: int(p.split('_')[1])):
        rows = [r for r in records if r['patient_id'] == pid]
        outcomes.append({'patient_id': pid, 'genes': sorted({r['gene'] for r in rows}),
                         'tmn_source_flag': int(pid in positive),
                         'clinical_records': clinical[pid],
                         'event_free_followup_available': False,
                         'label_interpretation': 'source_tMN_positive' if pid in positive else 'no_tMN_record_in_designated_source_not_confirmed_negative'})
    (DATA / 'outcome_evidence.json').write_text(json.dumps(outcomes, indent=2))
    audit = {'source_files': {p.name: {'sha256': hashlib.sha256(p.read_bytes()).hexdigest()} for p in (f7, f9)},
             'sheets_read': ['MOESM9:Gene_Level_Expansion_Data', 'MOESM7:tMN_VAF', 'MOESM7:Auto_CAR_tMN'],
             'longitudinal_observations_used': False,
             'total_variants': len(records), 'total_patients': len(ids),
             'source_positive_patients': len(ids & positive),
             'flags': dict(Counter(flag for r in records for flag in r['quality_flags'].split(';') if flag)),
             'genes': {}}
    for gene in GENES:
        rr = [r for r in records if r['gene'] == gene]
        eligible = [r for r in rr if 0 < r['first_vaf_percent'] < 50 and 0 < r['last_vaf_percent'] < 50]
        audit['genes'][gene] = {'variants': len(rr), 'patients': len({r['patient_id'] for r in rr}),
            'source_positive_patients': len({r['patient_id'] for r in rr if r['tmn_source_flag']}),
            'positive_diploid_endpoint_variants': len(eligible),
            'patient_gene_groups_with_positive_diploid_endpoints': len({r['patient_id'] for r in eligible})}
    (DATA / 'audit.json').write_text(json.dumps(audit, indent=2))
    print(json.dumps(audit, indent=2))

if __name__ == '__main__':
    main()
