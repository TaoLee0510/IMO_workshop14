"""Audit sources and prepare a patient-level, compartment-specific index-marker cohort."""
from __future__ import annotations
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
import openpyxl

ROOT = Path(__file__).resolve().parent
CONFIG = json.loads((ROOT / 'config.json').read_text())
DATA = ROOT / 'data'


def write_csv(path, rows, fields=None):
    fields = fields or list(rows[0])
    with path.open('w', newline='') as f:
        out = csv.DictWriter(f, fieldnames=fields)
        out.writeheader()
        out.writerows(rows)


def valid_number(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool)


def main():
    DATA.mkdir(exist_ok=True)
    w9path = ROOT.parent / '41588_2026_2795_MOESM9_ESM.xlsx'
    w7path = ROOT.parent / '41588_2026_2795_MOESM7_ESM.xlsx'
    old = list(csv.DictReader((ROOT.parent / 'moran_tmn_gt2_study/data/all_endpoints.csv').open()))
    genes = set(CONFIG['genes'])
    candidates = [r for r in old if r['gene'] in genes]
    loci = defaultdict(list)
    for r in candidates:
        chrom, pos, alleles = r['variant_id'].split(':')
        a, b = alleles.split('>')
        loci[(r['patient_id'], r['gene'], chrom, int(pos), tuple(sorted((a, b))))].append(r)
    w9 = openpyxl.load_workbook(w9path, read_only=True, data_only=True)
    metadata = defaultdict(list)
    for excelrow, row in enumerate(w9['Sample_meta_Data'].values, 1):
        if excelrow == 1:
            continue
        metadata[(row[0], row[1])].append(dict(row=excelrow, sample_type=row[2], day=row[5], depth=row[3]))
    observations = []
    scanned = 0
    orientation = Counter()
    for excelrow, row in enumerate(w9['Filtered_Mutation_Data'].values, 1):
        if excelrow == 1:
            continue
        scanned += 1
        if row[6] not in genes or not valid_number(row[3]):
            continue
        key = (row[0], row[6], row[2], int(row[3]), tuple(sorted((str(row[4]), str(row[5])))))
        matches = loci.get(key, [])
        for candidate in matches:
            flags = []
            if len(matches) != 1:
                flags.append('ambiguous_endpoint_locus_match')
            md = metadata.get((row[0], row[1]), [])
            if len(md) != 1:
                flags.append('missing_or_duplicate_sample_metadata')
            if len(md) == 1 and md[0]['day'] != row[10]:
                flags.append('sample_metadata_day_mismatch')
            if len(md) == 1 and md[0]['sample_type'] != row[11]:
                flags.append('sample_metadata_type_mismatch')
            source_alleles = candidate['variant_id'].split(':')[-1]
            direct = source_alleles == f'{row[4]}>{row[5]}'
            orient = 'TUMOR_to_REF' if direct else 'REF_to_TUMOR'
            orientation[orient] += 1
            v = row[8] if valid_number(row[8]) else None
            if v is None or not 0 < v < 50:
                flags.append('outside_positive_diploid_mapping')
            observations.append(dict(patient_id=row[0], gene=row[6], variant_id=candidate['variant_id'],
                sample_id=row[1], day=row[10], sample_type=row[11], vaf_percent=v,
                source_sheet='MOESM9:Filtered_Mutation_Data', source_excel_row=excelrow,
                endpoint_source_excel_row=candidate['source_excel_row'], allele_orientation=orient,
                metadata_excel_rows=';'.join(str(m['row']) for m in md), quality_flags=';'.join(flags)))
    w9.close()
    w7 = openpyxl.load_workbook(w7path, read_only=True, data_only=True)
    clinical = defaultdict(list)
    for excelrow, row in enumerate(w7['Auto_CAR_tMN'].values, 1):
        if excelrow > 1 and row[1] == 'CAR':
            clinical[row[0]].append(dict(source_excel_row=excelrow, time_to_tmn=row[2], multiple_CAR=row[11]))
    positives = set()
    for excelrow, row in enumerate(w7['tMN_VAF'].values, 1):
        if excelrow > 1 and row[-1] == 'tMN':
            positives.add(row[0])
    w7.close()
    pids = sorted({r['patient_id'] for r in old}, key=lambda x: int(x.split('_')[1]))
    outcomes = []
    for pid in pids:
        rows = clinical[pid]
        days = [r['time_to_tmn'] for r in rows if valid_number(r['time_to_tmn'])]
        flags = []
        if len(rows) > 1:
            flags.append('multiple_CAR_time_origins')
        if rows and pid not in positives:
            flags.append('clinical_event_without_designated_source_flag')
        if pid in positives and not days:
            flags.append('source_positive_missing_clinical_time')
        if days and min(days) > 1460:
            flags.append('listed_event_after_four_years')
        if any(str(r['multiple_CAR']) not in ('0', '0.0', 'No', 'no') for r in rows):
            flags.append('multiple_CAR_field_requires_review')
        outcomes.append(dict(patient_id=pid, source_positive=int(pid in positives),
            clinical_excel_rows=';'.join(str(r['source_excel_row']) for r in rows),
            candidate_event_days=';'.join(str(d) for d in days),
            clinical_status='unresolved' if rows or pid in positives else 'unknown',
            verified_event_day='', verified_censor_day='', verified_death_day='',
            index_infusion_verified=False, first_event_verified=False,
            ascertainment_verified=False, clinical_likelihood_ready=False,
            review_notes=';'.join(flags) or 'no_verified_event_free_followup'))
    # Duplicate sample/variant calls cannot be silently averaged or counted twice.
    duplicates = Counter((r['patient_id'], r['variant_id'], r['sample_id']) for r in observations)
    for r in observations:
        if duplicates[(r['patient_id'], r['variant_id'], r['sample_id'])] > 1:
            r['quality_flags'] = ';'.join(filter(None, [r['quality_flags'], 'duplicate_variant_sample_call']))
    bymarker = defaultdict(list)
    for r in observations:
        if r['sample_type'] in CONFIG['observation']['primary_sample_types']:
            bymarker[(r['patient_id'], r['gene'], r['variant_id'], r['sample_type'])].append(r)
    options = defaultdict(list)
    for key, rr in bymarker.items():
        baseline = [r for r in rr if valid_number(r['day']) and r['day'] <= 0
                    and valid_number(r['vaf_percent']) and 0 < r['vaf_percent'] < 50
                    and 'duplicate_variant_sample_call' not in r['quality_flags']
                    and 'sample_metadata_day_mismatch' not in r['quality_flags']]
        if not baseline:
            continue
        latest = max(r['day'] for r in baseline)
        atlatest = [r for r in baseline if r['day'] == latest]
        if len(atlatest) != 1:
            continue
        base = atlatest[0]
        options[key[0]].append((key, base, rr))
    selected, cohort = [], []
    for pid in pids:
        clinical_days = [r['time_to_tmn'] for r in clinical[pid] if valid_number(r['time_to_tmn'])]
        if len(clinical[pid]) > 1 or (pid in positives and not clinical_days):
            reason = 'unresolved_multiple_CAR_origins' if len(clinical[pid]) > 1 else 'known_positive_missing_diagnosis_time'
            cohort.append(dict(patient_id=pid, included=False, reason=reason,
                               gene='', variant_id='', sample_type='', baseline_day='', baseline_vaf_percent='', followup_observations=0))
            continue
        if not options[pid]:
            cohort.append(dict(patient_id=pid, included=False, reason='no_unique_positive_pre_CAR_blood_cell_baseline',
                               gene='', variant_id='', sample_type='', baseline_day='', baseline_vaf_percent='', followup_observations=0))
            continue
        # Baseline-only choice; follow-up and outcome do not influence marker ranking.
        key, base, rr = sorted(options[pid], key=lambda x: (-x[1]['vaf_percent'], x[0][1], x[0][2], x[0][3]))[0]
        eventdays = [r['time_to_tmn'] for r in clinical[pid] if valid_number(r['time_to_tmn'])]
        # Earliest candidate used only as a conservative VAF exclusion boundary, not an adjudicated event label.
        stop = min([1460.000001] + eventdays)
        usable = [r for r in rr if valid_number(r['day']) and base['day'] <= r['day'] < stop
                  and valid_number(r['vaf_percent']) and 0 < r['vaf_percent'] < 50
                  and 'duplicate_variant_sample_call' not in r['quality_flags']
                  and 'sample_metadata_day_mismatch' not in r['quality_flags']
                  and 'ambiguous_endpoint_locus_match' not in r['quality_flags']]
        daycounts = Counter(r['day'] for r in usable)
        usable = [r for r in usable if daycounts[r['day']] == 1]
        later = [r for r in usable if r['day'] > base['day']]
        include = len(later) > 0 and any(r['day'] == base['day'] for r in usable)
        cohort.append(dict(patient_id=pid, included=include, reason='eligible_VAF_diagnostic' if include else 'selected_baseline_marker_has_no_eligible_followup',
            gene=key[1], variant_id=key[2], sample_type=key[3], baseline_day=base['day'], baseline_vaf_percent=base['vaf_percent'], followup_observations=len(later)))
        if include:
            for r in sorted(usable, key=lambda x: x['day']):
                selected.append(dict(**r, role='initial_condition' if r['day'] == base['day'] else 'likelihood_observation'))
    write_csv(DATA / 'longitudinal_matched.csv', observations)
    write_csv(DATA / 'clinical_evidence_inventory.csv', outcomes)
    # Preserve human-adjudicated fields across source re-extraction.
    if not (DATA / 'clinical_reconciliation.csv').exists():
        write_csv(DATA / 'clinical_reconciliation.csv', outcomes)
    write_csv(DATA / 'index_marker_inventory.csv', cohort)
    if selected:
        write_csv(DATA / 'analysis_longitudinal.csv', selected)
    audit = dict(source_files={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in (w7path, w9path)},
        scanned_mutation_rows=scanned, matched_observation_rows=len(observations),
        matched_sample_types=dict(Counter(r['sample_type'] for r in observations)),
        allele_orientation=dict(orientation), quality_flags=dict(Counter(f for r in observations for f in r['quality_flags'].split(';') if f)),
        source_patients=len(pids), selected_index_patients=sum(r['included'] for r in cohort),
        selected_observations=len(selected), selected_likelihood_observations=sum(r['role']=='likelihood_observation' for r in selected),
        selected_sample_types=dict(Counter(r['sample_type'] for r in selected)),
        cohort_exclusions=dict(Counter(r['reason'] for r in cohort if not r['included'])),
        source_positive_patients=len(set(pids)&positives), clinical_ready_patients=0,
        clinically_verified_negative_patients=0, full_clinical_joint_fit_ready=False,
        clinical_missing_fields=['verified index infusion', 'first diagnosis adjudication', 'event-free censor dates', 'death dates', 'ascertainment review'],
        notes=['Unknown clinical outcomes are not negative.', 'Baseline-only marker selection precedes follow-up eligibility.',
               'PBL/PBMC not pooled within a patient; cfDNA/tumor excluded from primary blood-cell mapping.',
               'Patients with multiple CAR origins or source-positive outcomes without a diagnosis time are excluded from the initial primary VAF diagnostic.',
               'Positive-call table is not a complete assay detection matrix; absent rows are not zero.',
               'Earliest unadjudicated event time is used only to conservatively exclude potentially postdiagnosis VAF.'])
    (DATA / 'audit.json').write_text(json.dumps(audit, indent=2))
    print(json.dumps(audit, indent=2))


if __name__ == '__main__':
    main()
