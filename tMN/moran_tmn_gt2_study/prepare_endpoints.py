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
CONFIG = json.loads((ROOT / 'approved_scenarios.json').read_text())
GENES = tuple(CONFIG['genes'])

def eligible(record):
    return (all(isinstance(record[k], (int, float)) for k in ('first_day', 'last_day', 'first_vaf_percent', 'last_vaf_percent')) and
            record['first_day'] < record['last_day'] and
            0 < record['first_vaf_percent'] < 50 and
            0 < record['last_vaf_percent'] < 50)

def write_csv(path, records):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)

def main():
    f7 = SOURCE / '41588_2026_2795_MOESM7_ESM.xlsx'
    f9 = SOURCE / '41588_2026_2795_MOESM9_ESM.xlsx'
    w7 = openpyxl.load_workbook(f7, read_only=True, data_only=True)
    vaf_rows = list(w7['tMN_VAF'].values)
    assert vaf_rows[0][0] == 'Paper_ID' and vaf_rows[0][-1] == 'Timpoint'
    positive = {r[0] for r in vaf_rows[1:] if r[-1] == 'tMN'}
    positive_rows = defaultdict(list)
    for excel_row, row in enumerate(vaf_rows[1:], 2):
        if row[-1] == 'tMN':
            positive_rows[row[0]].append(excel_row)
    clinical = defaultdict(list)
    for excel_row, row in enumerate(w7['Auto_CAR_tMN'].values, 1):
        if excel_row > 1:
            clinical[row[0]].append({'source_excel_row': excel_row, 'cohort': row[1], 'time_cell_therapy_to_tMN': row[2]})
    w7.close()
    w9 = openpyxl.load_workbook(f9, read_only=True, data_only=True)
    # Do not read Filtered_Mutation_Data or any intermediate VAF measurements.
    records = []
    for excel_row, row in enumerate(w9['Gene_Level_Expansion_Data'].values, 1):
        if excel_row == 1 or not any(value is not None for value in row):
            continue
        assert row[1] == row[8]
        assert row[2] < row[3]
        flags = []
        first_vaf = row[4] if isinstance(row[4], (int, float)) else None
        last_vaf = row[5] if isinstance(row[5], (int, float)) else None
        if first_vaf is None:
            flags.append('missing_baseline_VAF')
        if last_vaf is None:
            flags.append('missing_last_VAF')
        if row[4] == 0:
            flags.append('baseline_zero_requires_detection_model')
        if row[5] == 0:
            flags.append('last_zero_requires_detection_model')
        if any(v is not None and v >= 50 for v in (first_vaf, last_vaf)):
            flags.append('diploid_heterozygous_mapping_not_valid_for_fit')
        if first_vaf is not None and last_vaf is not None and last_vaf < first_vaf:
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
            'first_vaf_percent': first_vaf, 'last_vaf_percent': last_vaf,
            'first_vaf_source_value': row[4], 'last_vaf_source_value': row[5],
            'source_excel_row': excel_row,
            'tmn_source_flag': int(row[0] in positive),
            'tmn_status_interpretation': 'source_tMN_positive' if row[0] in positive else 'no_tMN_record_not_confirmed_clinical_negative',
            'tmn_source_excel_rows': ';'.join(map(str, positive_rows[row[0]])),
            'quality_flags': ';'.join(flags),
        })
    w9.close()
    assert len(records) == 171
    assert len({(r['patient_id'], r['gene'], r['variant_id']) for r in records}) == len(records)
    gene_inventory = []
    for gene in sorted({r['gene'] for r in records}):
        rr = [r for r in records if r['gene'] == gene]
        usable = [r for r in rr if eligible(r)]
        count = len({r['patient_id'] for r in usable})
        gene_inventory.append({'gene': gene, 'all_variant_records': len(rr),
            'all_patient_gene_groups': len({r['patient_id'] for r in rr}),
            'eligible_variant_records': len(usable), 'eligible_patient_gene_groups': count,
            'retained_by_gt2_filter': count > 2,
            'status': 'retained' if count > 2 else 'excluded_eligible_patient_count_not_above_two',
            'source_tMN_positive_patients': len({r['patient_id'] for r in rr if r['tmn_source_flag']})})
    retained = {r['gene'] for r in gene_inventory if r['retained_by_gt2_filter']}
    assert retained == set(GENES), 'Gene inventory differs from the confirmed six-gene scope.'
    for record in records:
        record['gene_retained'] = record['gene'] in retained
        record['eligible_endpoint_pair'] = eligible(record)
        record['analysis_included'] = record['gene_retained'] and record['eligible_endpoint_pair']
    selected = [r for r in records if r['gene_retained']]
    assert len(selected) == 155
    assert sum(r['analysis_included'] for r in records) == 110
    assert len({(r['patient_id'], r['gene']) for r in records if r['analysis_included']}) == 83
    write_csv(DATA / 'all_endpoints.csv', records)
    write_csv(DATA / 'endpoints.csv', selected)
    write_csv(DATA / 'gene_eligibility.csv', gene_inventory)
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
             'gene_filter': CONFIG['gene_filter'],
             'total_variants': len(records), 'total_patients': len(ids),
             'all_genes': len(gene_inventory), 'all_patient_gene_groups': len({(r['patient_id'], r['gene']) for r in records}),
             'retained_genes': list(GENES), 'selected_variants': len(selected),
             'selected_patient_gene_groups': len({(r['patient_id'], r['gene']) for r in selected}),
             'eligible_selected_variants': sum(r['analysis_included'] for r in records),
             'eligible_selected_patient_gene_groups': sum(r['eligible_patient_gene_groups'] for r in gene_inventory if r['retained_by_gt2_filter']),
             'source_positive_patients': len(ids & positive),
             'flags': dict(Counter(flag for r in records for flag in r['quality_flags'].split(';') if flag)),
             'genes': {}}
    for gene in sorted({r['gene'] for r in records}):
        rr = [r for r in records if r['gene'] == gene]
        eligible_rows = [r for r in rr if eligible(r)]
        audit['genes'][gene] = {'variants': len(rr), 'patients': len({r['patient_id'] for r in rr}),
            'source_positive_patients': len({r['patient_id'] for r in rr if r['tmn_source_flag']}),
            'positive_diploid_endpoint_variants': len(eligible_rows),
            'patient_gene_groups_with_positive_diploid_endpoints': len({r['patient_id'] for r in eligible_rows}),
            'retained_by_gt2_filter': gene in retained}
    (DATA / 'audit.json').write_text(json.dumps(audit, indent=2))
    print(json.dumps(audit, indent=2))

if __name__ == '__main__':
    main()
