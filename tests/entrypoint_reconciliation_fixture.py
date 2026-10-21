"""Saved-file adapters for the UNCHANGED main viSUM workflow.

Only VIHARMONY remains a process. All other imported modules are inert channel
adapters: no tool binary, preparation, standardizer, or preflight is executed.
This checks main-workflow joins/reference wiring, NOT upstream computation.
"""
from __future__ import annotations
from collections import Counter
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import shutil

BUNDLE = 'assets/ictv_reconciliation/MSL41.v1.e83bc2c73a74.json.gz'
SAMPLES = ('BMARK_DNA_CELLULAR', 'BMARK_RNA_VIRAL')
TOOLS = ('genomad', 'virsorter2', 'cenotetaker3', 'deep6', 'deepmicroclass2',
         'virbot', 'gianthunter', 'vicat', 'checkv', 'tesorter', 'vitap', 'vcontact3')
BASE_MANIFEST = '948cc71190f8d079d169586c14da6607e260254e8398868ad00fb7902b2ae5ef'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding='utf-8', newline='\n')


def verify_original(project):
    project = project.resolve()
    manifest_path = project / 'PACKAGE_MANIFEST.json'
    if sha(manifest_path) != BASE_MANIFEST:
        raise ValueError('Not the original verified integrated package')
    manifest = json.loads(manifest_path.read_text())
    for name, digest in manifest['files'].items():
        path = (project / name).resolve()
        if project not in path.parents or sha(path) != digest:
            raise ValueError('Original package changed: ' + name)
    return manifest


def specs(original, project):
    """Use exactly the earlier native gate's saved DNA/RNA inputs.

    The old snapshots lack a separate viCAT context table. A clearly identified
    header-only neutral fixture occupies that current channel slot; it contains
    no prediction rows. DNA loci are likewise empty and never fed to viHARMONY.
    """
    result = []
    for sample in SAMPLES:
        source = original / 'saved/inputs' / (sample + '_results')
        fasta = original / 'saved/fasta' / (sample + '_results')
        published = original / 'saved/published' / sample
        old = json.loads((published / 'updated' / (sample + '.harmonizer_manifest.json')).read_text())
        evidence = {}
        for tool in TOOLS:
            path = published / (sample + '.virsorter2_evidence.tsv') if tool == 'virsorter2' else source / tool / (sample + '.' + tool + '_evidence.tsv')
            if path.is_file():
                evidence[tool] = str(path)
        expected_tools = set(TOOLS) - ({'deep6', 'virbot'} if old['input_type'] == 'dna' else {'deepmicroclass2', 'gianthunter'})
        if set(evidence) != expected_tools:
            raise ValueError('Missing/extra saved tool evidence: ' + sample)
        context = project / 'fixture-data' / (sample + '.vicat_context.tsv')
        write(context, 'sample_id\tsequence_id\tparent_sequence_id\ttool\tclassification\n')
        loci = source / 'vicat' / (sample + '.vicat_orf_evidence.tsv')
        if old['input_type'] == 'dna':
            loci = project / 'fixture-data' / loci.name
            write(loci, 'sample_id\tsequence_id\tcoordinates\tlocus_classification\tprotein_length\tbest_bitscore\tbest_query_coverage\tbest_subject_coverage\n')
        paths = dict(normalized=fasta / 'prep' / (sample + '.normalized.fasta'),
            refined=fasta / 'refinement' / (sample + '.refined_candidates.fasta'),
            header=source / 'prep' / (sample + '.header_map.tsv'),
            gate=source / 'discovery_gate' / (sample + '.discovery_gate.tsv'),
            regions=source / 'refinement' / (sample + '.provirus_region_map.tsv'),
            groups=source / 'vcontact3' / (sample + '.vcontact3_group_membership.tsv'),
            projected=source / 'vicat' / (sample + '.vicat_refined_evidence.tsv'),
            context=context, loci=loci)
        if not all(p.is_file() for p in paths.values()):
            raise ValueError('Missing saved input: ' + sample)
        result.append(dict(sample=sample, type=old['input_type'], evidence=evidence,
            **{key: str(path) for key, path in paths.items()},
            minimum_length=old['policy']['vcontact3_min_taxonomy_length'],
            pair_floor=old['policy']['rna_pair_homology_floor']))
    return result


def adapter(name, inputs, outputs, special=None):
    take = '\n'.join('    arg' + str(i) for i in range(inputs))
    main, emit = [], []
    for label, expression in outputs.items():
        if special == 'preflight':
            channel = 'arg0.map { row -> row[0] }'
        elif special == 'reference':
            channel = "arg0.map { ignored -> def p = file(params.entry_specs[0].normalized, checkIfExists: true); tuple(p, p) }"
        else:
            channel = ("arg0.map { row ->\n"
                "        def s = params.entry_specs.find { it.sample == row[0] }\n"
                "        if( s == null ) error 'Unexpected sample in saved-input adapter'\n"
                "        def p = file(s.normalized, checkIfExists: true)\n"
                "        " + expression + '\n    }')
        main.append('    out_' + label + ' = ' + channel)
        emit.append('    ' + label + ' = out_' + label)
    return ('// TEST ONLY: saved-file channel adapter; no external command.\nworkflow ' + name + ' {\n'
            '    take:\n' + take + '\n    main:\n' + '\n'.join(main) + '\n    emit:\n' + '\n'.join(emit) + '\n}\n')


def adapters():
    records = {}
    def add(name, inputs, outputs, special=None):
        records[name] = adapter(name, inputs, outputs, special)
    def raw(name, inputs, arity):
        add(name, inputs, dict(results='tuple(row[0], row[1], ' + ', '.join(['p'] * (arity - 2)) + ')'))
    def evidence(tool, inputs=1):
        add('STANDARDIZE_' + tool.upper(), inputs, dict(evidence="tuple(row[0], '" + tool + "', file(s.evidence['" + tool + "'], checkIfExists: true))"))
    add('VALIDATE_PREPARED_ASSETS', 1, dict(samples=''), 'preflight')
    add('NORMALIZE_FASTA', 1, dict(normalized_records="tuple(row[0], row[1], p, file(s.header, checkIfExists: true))"))
    for tool, arity in [('genomad', 9), ('virsorter2', 7), ('cenotetaker3', 10),
                        ('deep6', 5), ('deepmicroclass2', 5), ('virbot', 6), ('gianthunter', 6),
                        ('checkv', 9)]:
        raw('RUN_' + tool.upper(), 2, arity)
        evidence(tool, 2 if tool in ('virbot', 'gianthunter') else 1)
    add('PREDICT_VICAT_ORFS', 1, dict(orfs="tuple(row[0], row[1], p, p, file(s.header, checkIfExists: true))"))
    raw('RUN_VICAT_DIAMOND', 2, 5)
    raw('RUN_VICAT_NONVIRAL_DIAMOND', 2, 5)
    add('PREPARE_VICAT_REFERENCE_SUBSET', 2, dict(subset=''), 'reference')
    raw('PREPARE_VICAT_HITS', 3, 6)
    add('STANDARDIZE_VICAT', 2, dict(
        evidence="tuple(row[0], 'vicat', file(s.evidence.vicat, checkIfExists: true))",
        loci="tuple(row[0], file(s.loci, checkIfExists: true))",
        context="tuple(row[0], 'vicat_context', file(s.context, checkIfExists: true))",
        provirus="tuple(row[0], 'vicat', file(s.evidence.vicat, checkIfExists: true))"))
    add('PROJECT_VICAT_REFINED', 1, dict(
        evidence="tuple(row[0], 'vicat', file(s.projected, checkIfExists: true))",
        projection="tuple(row[0], file(s.projected, checkIfExists: true))"))
    add('DISCOVERY_GATE', 1, dict(candidates="tuple(row[0], row[1], p, file(s.gate, checkIfExists: true), p)"))
    add('REFINE_PROVIRAL_REGIONS', 2, dict(refined="tuple(row[0], row[1], file(s.refined, checkIfExists: true), file(s.regions, checkIfExists: true), p, p)"))
    add('RUN_TESORTER', 1, dict(results="tuple(row[0], row[1], file(s.regions, checkIfExists: true), p, p, p, p, p, p)",
                              completed='tuple(row[0], row[1])'))
    evidence('tesorter')
    raw('RUN_VITAP', 3, 8)
    evidence('vitap', 2)
    raw('RUN_VCONTACT3', 3, 7)
    add('STANDARDIZE_VCONTACT3', 1, dict(
        evidence="tuple(row[0], 'vcontact3', file(s.evidence.vcontact3, checkIfExists: true))",
        groups="tuple(row[0], file(s.groups, checkIfExists: true))"))
    return records


def imports(text):
    return re.findall(r"include \{ (\w+) \} from '([^']+)'", text)


def tuple_arity(text):
    """Count outer tuple fields, ignoring nested file(...) argument commas."""
    start = text.index('tuple(') + len('tuple(')
    depth, quote, fields = 1, None, 1
    for char in text[start:]:
        if quote:
            if char == quote:
                quote = None
        elif char in "\"'":
            quote = char
        elif char == '(':
            depth += 1
        elif char == ')':
            depth -= 1
            if depth == 0:
                return fields
        elif char == ',' and depth == 1:
            fields += 1
    raise ValueError('Unterminated tuple expression')


def interface_contracts(source_root):
    """Check every adapter's input count and emitted tuple arity against real modules."""
    result = {}
    include_map = dict(imports((source_root / 'visum_nextflow.nf').read_text()))
    for name, generated in adapters().items():
        module = (source_root / (include_map[name].removeprefix('./') + '.nf')).read_text()
        inputs = module.split('    input:', 1)[1].split('    output:', 1)[0]
        original_inputs = len(re.findall(r'^    (?:tuple|val|path)\b', inputs, flags=re.M))
        take = generated.split('    take:', 1)[1].split('    main:', 1)[0]
        if len(re.findall(r'^    arg\d+$', take, flags=re.M)) != original_inputs:
            raise AssertionError('Adapter input count differs from real module: ' + name)
        output = module.split('    output:', 1)[1].split('    script:', 1)[0]
        declared = {}
        for match in re.finditer(r'^    ((?:tuple|val|path)\b.*?emit:\s*(\w+))', output, flags=re.M | re.S):
            text, label = match.group(1), match.group(2)
            declared[label] = len(re.findall(r'\b(?:val|path)\(', text)) if text.startswith('tuple') else 1
        emitted = {}
        for match in re.finditer(r'^    out_(\w+) = (.*?)(?=\n    out_|\n    emit:)', generated, flags=re.M | re.S):
            label, body = match.group(1), match.group(2)
            arity = tuple_arity(body) if 'tuple(' in body else 1
            if declared.get(label) != arity:
                raise AssertionError('Adapter output differs from real module: ' + name + '/' + label)
            emitted[label] = arity
        result[name] = dict(inputs=original_inputs, used_outputs=emitted)
    return result


def build(original, output):
    """Generate ONLY inside a new scratch project; main/scientific files copied verbatim."""
    project = output / 'project'
    if project.exists():
        raise ValueError('Fixture project already exists; use a new output directory')
    project.mkdir(parents=True)
    for name in ('visum_nextflow.nf', 'visum.config', 'nextflow.config',
                 'modules/local/viharmony.nf', 'bin/run_viharmony.py',
                 'bin/evidence_schema.py', 'bin/ictv_reconciliation.py',
                 'bin/standardize_vicat.py', 'assets/ICTV_VMR_MSL41.csv', BUNDLE, BUNDLE + '.manifest.json'):
        destination = project / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(original / name, destination)
        if sha(original / name) != sha(destination):
            raise AssertionError('Copy not byte-identical: ' + name)
    write(project / 'bin/refine_proviral_regions.py', '# TEST PLACEHOLDER: never executed; refinement supplied by saved-file adapter.\n')
    write(project / 'workflows/database_setup.nf', "workflow DATABASE_SETUP {\n    main:\n    error 'STOP: Database setup is forbidden in the entry-point fixture'\n}\n")
    generated = adapters()
    included = imports((project / 'visum_nextflow.nf').read_text())
    if {name for name, path in included} != set(generated) | {'VIHARMONY', 'DATABASE_SETUP'}:
        raise ValueError('Main imports changed; review all saved-file adapters')
    for name, path in included:
        if name in generated:
            write(project / (path.removeprefix('./') + '.nf'), generated[name])
    selected = specs(original, project)
    (project / 'data').mkdir(exist_ok=True)
    # Absolute FASTA paths are permitted in the actual sample-sheet interface.
    sheet = project / 'data/samples.csv'
    with sheet.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['prefix', 'type', 'fasta'])
        writer.writerows((s['sample'], s['type'], s['normalized']) for s in selected)
    if len({s['minimum_length'] for s in selected}) != 1:
        raise ValueError('Selected panels require different taxonomy-length settings')
    payload = dict(entry_specs=selected, prefix_many=str(sheet), indir=str(project / 'data'),
        outdir=str(output / 'published'), setup=False, harmonizer_cpus=8, max_cpus=8,
        harmonizer_max_rank='species',  # compare to the frozen pre-cap reconciliation results
        harmonizer_memory='32 GB', max_memory='64 GB', harmonizer_time=None,
        harmonizer_audit='full', vcontact3_min_taxonomy_length=selected[0]['minimum_length'],
        rna_pair_homology_floor=next(s['pair_floor'] for s in selected if s['type'] == 'rna'),
        show_channel_messages=False)
    payload.update({'run_' + tool: True for tool in TOOLS})
    # No database is read. Preflight is explicitly a fixture adapter, not a test of readiness.
    (project / 'databases').mkdir()
    write(output / 'params.json', json.dumps(payload, indent=2) + '\n')
    protected = {name: sha(project / name) for name in (
        'visum_nextflow.nf', 'visum.config', 'nextflow.config', 'modules/local/viharmony.nf',
        'bin/run_viharmony.py', 'bin/evidence_schema.py', 'bin/ictv_reconciliation.py')}
    write(output / 'fixture_inventory.json', json.dumps(dict(
        scope='Actual main workflow + actual VIHARMONY; every other module is an inert saved-file adapter',
        original=str(original), protected_files=protected, adapted_modules=sorted(generated),
        neutral_fixture_slots=['viCAT context: no rows; DNA loci: no rows, not sent to harmonizer'],
        samples=selected), indent=2) + '\n')
    return project, selected, protected


def output_manifest(folder, sample):
    result = json.loads((folder / (sample + '.harmonizer_manifest.json')).read_text())
    for item in result['outputs'].values():
        if Path(item['file']).name != item['file'] or sha(folder / item['file']) != item['sha256']:
            raise AssertionError('Invalid output manifest/checksum: ' + str(folder))
    return result


def canonical_audit(path):
    data = gzip.decompress(path.read_bytes()).decode('utf-8-sig')
    reader = csv.DictReader(io.StringIO(data), delimiter='\t')
    fields = reader.fieldnames
    if len(fields) != len(set(fields)):
        raise AssertionError('Duplicate audit columns')
    return set(fields), Counter(json.dumps(row, sort_keys=True) for row in reader)


def compare_outputs(expected, actual, sample):
    """Exact scientific outputs; audit table ordering is not a biological difference.

    Real main workflow channel merges can change evidence arrival order. Only
    *_audit.tsv.gz row/column order and gzip timestamps are normalized; all
    columns, row multiplicities, values and other output bytes must match.
    """
    left, right = output_manifest(expected, sample), output_manifest(actual, sample)
    if left['outputs'].keys() != right['outputs'].keys():
        raise AssertionError('Output sets differ')
    for label, item in left['outputs'].items():
        other = right['outputs'][label]
        if item['file'] != other['file']:
            raise AssertionError('Output filename changed')
        a, b = expected / item['file'], actual / other['file']
        if a.name.endswith('_audit.tsv.gz'):
            same = canonical_audit(a) == canonical_audit(b)
        else:
            same = a.read_bytes() == b.read_bytes()
        if not same:
            raise AssertionError('Scientific output differs: ' + sample + '/' + label)
        item.pop('sha256')
        other.pop('sha256')
    if left != right:
        raise AssertionError('Policy/reference/output metadata changed: ' + sample)
    return dict(sample=sample, outputs_verified=len(left['outputs']), scientific_outputs_identical=True,
        audit_order_only_normalized=True)


def expected_evidence(spec):
    paths = list(spec['evidence'].values()) + [spec['projected'], spec['context']]
    if spec['type'] == 'rna':
        paths.append(spec['loci'])
    return paths


def verify_command(script, spec, reference):
    """Inspect the REAL task command: every saved evidence file and reference must arrive."""
    import shlex
    text = script.read_text()
    # Original module uses shell continuation lines and one python3 command.
    commands = [line for line in text.replace('\\\n', ' ').splitlines() if line.strip().startswith('python3 ')]
    if len(commands) != 1:
        raise AssertionError('Expected exactly one real harmonizer invocation')
    argv = shlex.split(commands[0])
    def value(flag):
        return argv[argv.index(flag) + 1]
    if value('--sample-id') != spec['sample'] or value('--input-type') != spec['type']:
        raise AssertionError('Main workflow mixed sample/type')
    for flag, key in [('normalized-fasta', 'normalized'), ('header-map', 'header'),
                      ('discovery-gate', 'gate'), ('refined-fasta', 'refined'), ('region-map', 'regions')]:
        if value('--' + flag) != Path(spec[key]).name:
            raise AssertionError('Main workflow passed incorrect input: ' + flag)
    start, stop = argv.index('--evidence') + 1, argv.index('--vcontact3-groups')
    if Counter(argv[start:stop]) != Counter(Path(p).name for p in expected_evidence(spec)):
        raise AssertionError('Main workflow dropped/duplicated/mixed evidence: ' + spec['sample'])
    if value('--vcontact3-groups') != Path(spec['groups']).name:
        raise AssertionError('Wrong vCONTACT3 groups')
    if value('--ictv-reconciliation') != reference.name or value('--workers') != '8':
        raise AssertionError('Wrong pinned bundle/worker setting')
    if (value('--ictv-msl') != 'ICTV_VMR_MSL41.csv' or value('--audit-mode') != 'full'
            or value('--vcontact3-min-taxonomy-length') != str(spec['minimum_length'])
            or ('--rna-pair-homology-floor' in argv) != (spec['type'] == 'rna' and spec['pair_floor'])
            or '--disable-rna-checkv-exception' in argv):
        raise AssertionError('Main workflow changed reference/audit/scientific settings')
    return dict(sample=spec['sample'], evidence_files_verified=stop-start, command_inputs_verified=True)
