"""Native configuration/resource probes only; never run the scientific workflow."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = 'resource-defaults-20261015-v2'
PROFILES = ('default', 'local', 'local_safe', 'local_full', 'vicat_build',
            'slurm', 'docker', 'singularity')
UNFINISHED_PROFILES = ('docker', 'singularity')
PROBE_CASES = ('defaults', 'parameters', 'site')
GIB = 1024 ** 3


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(value, encoding='utf-8', newline='\n')


def time_parameters(root=ROOT):
    text = (root / 'visum.config').read_text(encoding='utf-8')
    pairs = re.findall(r'(?m)^\s*(\w+_time)\s*=\s*([^\n]+)', text)
    if len(pairs) != 20 or len(dict(pairs)) != 20:
        raise ValueError('Expected exactly twenty distinct wall-time parameters')
    if any(value.strip() != 'null' for _, value in pairs):
        raise ValueError('A shipped time default is not null')
    return tuple(name for name, _ in pairs)


def verify_snapshot(root=ROOT):
    snapshot = json.loads((root / 'RESOURCE_TEST_SNAPSHOT.json').read_text())
    if snapshot['test_package'] != PACKAGE:
        raise ValueError('Wrong frozen test package')
    for name, expected in snapshot['files'].items():
        path = root / name
        if not path.resolve().is_relative_to(root.resolve()) or path.is_symlink():
            raise ValueError(f'Unsafe snapshot member: {name}')
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError(f'Snapshot file changed: {name}')
    print(f"Verified snapshot: {len(snapshot['files'])} files", flush=True)
    return snapshot


def probe_workflow(root=ROOT):
    """Copy time expressions, not executable analysis/install bodies."""
    names = {'RESOURCE_UTILITY': 'utility_time'}
    processes = []
    definitions = [('RESOURCE_UTILITY', '', "memory '64 MB'\n    cpus 1")]
    for path in sorted((root / 'modules/local').glob('*.nf')):
        source = path.read_text(encoding='utf-8')
        times = re.findall(r'(?m)^\s*time\s+([^\n]+)', source)
        if not times:
            continue
        if len(times) != 1:
            raise ValueError(f'Unexpected time declarations: {path.name}')
        match = re.fullmatch(r'(?:params\.(\w+)|\{\s*params\.(\w+)\s*\})', times[0].strip())
        process = re.search(r'(?m)^process\s+(\w+)\s*\{', source)
        if not match or not process:
            raise ValueError(f'Unsupported resource declaration: {path.name}')
        name = 'RESOURCE_' + process[1]
        if name in names:
            raise ValueError(f'Duplicate probe: {name}')
        names[name] = match[1] or match[2]
        resources = "memory '64 MB'\n    cpus 1"
        if process[1] == 'VIHARMONY':
            # Fail closed if the real module's tested allocations drift.
            expected = ('memory params.harmonizer_memory',
                        'cpus { Math.min(params.harmonizer_cpus as int, params.max_cpus as int) }')
            if not all(line in source for line in expected):
                raise ValueError('viHARMONY resource declaration changed')
            resources = '\n    '.join(expected)
        definitions.append((name, 'time ' + times[0].strip(), resources))
    if len(names) < 26 or set(names.values()) != set(time_parameters(root)):
        raise ValueError('Probe graph does not cover every time parameter')
    for name, directive, resources in definitions:
        processes.append('''process NAME {
    RESOURCES
    TIME
    publishDir params.probe_outdir, mode: 'copy'
    output:
    path 'NAME.json'
    script:
    def observation = groovy.json.JsonOutput.toJson([
        probe: 'NAME', cpus: task.cpus,
        memory_bytes: task.memory.toBytes(),
        time_ms: task.time == null ? null : task.time.toMillis()
    ])
    """
    printf '%s\\n' '${observation}' > NAME.json
    """
}
'''.replace('NAME', name).replace('RESOURCES', resources).replace('TIME', directive))
    graph = '\n'.join(processes) + '\nworkflow {\n' + ''.join(
        f'    {name}()\n' for name in names) + '}\n'
    return graph, names


def parse_config_output(output):
    """Decode one JSON object, retaining WARN lines without hiding other output."""
    remaining = output.lstrip()
    warnings = []
    while remaining.startswith('WARN:'):
        line, separator, rest = remaining.partition('\n')
        warnings.append(line.rstrip('\r'))
        remaining = rest.lstrip() if separator else ''
    # raw_decode rejects malformed/missing JSON; do not search for a later '{'
    # that might conceal an ERROR, traceback, or unrelated console output.
    config, end = json.JSONDecoder().raw_decode(remaining)
    if not isinstance(config, dict):
        raise ValueError('Nextflow configuration must be a JSON object')
    for line in remaining[end:].splitlines():
        if not line.strip():
            continue
        if not line.startswith('WARN:'):
            raise ValueError('Unexpected output after configuration JSON: ' + line)
        warnings.append(line.rstrip('\r'))
    return config, warnings


def validate_config(config, profile, parameters):
    values = config['params']
    missing = set(parameters) - set(values)
    if missing or any(values[name] is not None for name in parameters):
        raise AssertionError(f'{profile}: missing/non-null time defaults: {missing}')
    if values['harmonizer_memory'] != '32 GB' or values['harmonizer_cpus'] != 1:
        raise AssertionError(f'{profile}: incorrect harmonizer defaults')
    expected_memory = {'local_full': '128 GB', 'vicat_build': '400 GB'}.get(profile, '64 GB')
    if values['max_memory'] != expected_memory:
        raise AssertionError(f'{profile}: incorrect scheduling memory budget')
    # JSON may render the deferred utility closure differently across parsers.
    # Reject literal durations here; the probes evaluate the closure itself.
    directive = config.get('process', {}).get('time')
    if directive is not None and re.match(r'^\s*\d', str(directive)):
        raise AssertionError(f'{profile}: finite process.time: {directive}')
    if 'pull' in config.get('docker', {}):
        raise AssertionError(f'{profile}: unsupported docker.pull setting remains')


def validate_observations(directory, names, case):
    actual = {path.stem: json.loads(path.read_text()) for path in directory.glob('*.json')}
    if set(actual) != set(names):
        raise AssertionError(f'{case}: missing/extra resource observations')
    for name, observation in actual.items():
        expected_time = 3 * 3600 * 1000 if case == 'parameters' else None
        harmony = name == 'RESOURCE_VIHARMONY'
        if case == 'site':
            if harmony:
                expected_time = 7 * 3600 * 1000
            elif name == 'RESOURCE_UTILITY':
                expected_time = 5 * 3600 * 1000
        expected = dict(probe=name, cpus=2 if harmony and case != 'defaults' else 1,
                        memory_bytes=(32 if case == 'defaults' else 48) * GIB if harmony else 64 * 1024 ** 2,
                        time_ms=expected_time)
        if observation != expected:
            raise AssertionError(f'{case}/{name}: {observation} != {expected}')
    return actual


def run_command(command, directory, label):
    directory.mkdir(parents=True, exist_ok=True)
    write(directory / (label + '.command.json'), json.dumps(command, indent=2) + '\n')
    # Safety watchdog belongs to this tiny fixture, not the shipped processes.
    with (directory / (label + '.console.log')).open('w', encoding='utf-8') as output:
        process = subprocess.Popen(command, cwd=directory, stdout=output,
                                   stderr=subprocess.STDOUT, start_new_session=True)
        try:
            status = process.wait(timeout=180)
        except subprocess.TimeoutExpired:
            import signal
            os.killpg(process.pid, signal.SIGTERM)  # This newly spawned fixture group only.
            try:
                process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait(timeout=10)
            raise RuntimeError(f'Fixture watchdog expired: {directory}')
    if status:
        raise RuntimeError(f'Nextflow exited {status}: {directory / (label + ".console.log")}')
    return (directory / (label + '.console.log')).read_text(encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify-only', action='store_true')
    parser.add_argument('--run-dir', type=Path)
    args = parser.parse_args()
    verify_snapshot()
    if args.verify_only:
        return 0
    if sys.platform != 'linux' or args.run_dir is None:
        parser.error('Use the packaged Linux runner; no tests were run')
    run = args.run_dir.resolve()
    if not run.is_relative_to(ROOT.parent.resolve()) or not run.name.startswith('resource-defaults-run.'):
        parser.error('Run directory must be a fresh fixture beside this snapshot')
    report = dict(test_package=PACKAGE, status='FAIL', profiles={}, probes={}, warnings={},
                  selection=dict(profiles=list(UNFINISHED_PROFILES), probes=list(PROBE_CASES)),
                  previously_passed_profiles_not_rerun=[p for p in PROFILES if p not in UNFINISHED_PROFILES],
                  scope='resource probes only; no scientific workflow or database access')
    started = time.monotonic()
    try:
        parameters = time_parameters()
        executable = str(run / 'bin/nextflow')
        if not Path(executable).is_file():
            raise ValueError('Private Nextflow launcher missing')
        for profile in UNFINISHED_PROFILES:
            directory = run / 'artifacts/config' / profile
            command = [executable, '-log', str(directory / 'nextflow.log'),
                       '-C', str(ROOT / 'nextflow.config'), 'config', str(ROOT), '-o', 'json']
            if profile != 'default':
                command += ['-profile', profile]
            output = run_command(command, directory, 'config')
            config, warnings = parse_config_output(output)
            write(directory / 'config.parsed.json', json.dumps(config, indent=2) + '\n')
            write(directory / 'config.warnings.json', json.dumps(warnings, indent=2) + '\n')
            report['warnings'][profile] = warnings
            for warning in warnings:
                print(f'{profile}: {warning}', flush=True)
            validate_config(config, profile, parameters)
            report['profiles'][profile] = 'PASS'
            print(f'PASS config: {profile}', flush=True)
        graph, names = probe_workflow()
        for case in PROBE_CASES:
            directory = run / 'artifacts/probes' / case
            project = directory / 'project'
            project.mkdir(parents=True)
            write(project / 'resource_probes.nf', graph)
            private = directory / 'observations'
            # Disable dependencies and nonlocal execution for EVERY probe.
            # Keep shipped memory/time expressions for viHARMONY; other probes
            # use 64 MB because their memory sizing is outside this test's scope.
            config = '''conda.enabled = false
docker.enabled = false
singularity.enabled = false
dag.enabled = false
report.enabled = false
timeline.enabled = false
trace.enabled = false
process {
    withName: 'RESOURCE_.*' {
        executor = 'local'
        errorStrategy = 'terminate'
        maxRetries = 0
    }
}
'''
            if case == 'site':
                config += '''process {
    withName: 'RESOURCE_UTILITY' { time = '5h' }
    withName: 'RESOURCE_VIHARMONY' {
        time = '7h'
        memory = '48 GB'
        cpus = 2
    }
}
'''
            write(project / 'probe.config', config)
            values = dict(max_cpus=2, max_memory='64 GB', probe_outdir=str(private))
            if case == 'parameters':
                values.update({key: '3h' for key in parameters})
                values.update(harmonizer_memory='48 GB', harmonizer_cpus=8)
            write(project / 'params.json', json.dumps(values, indent=2) + '\n')
            command = [executable, '-log', str(directory / 'nextflow.log'),
                       '-C', str(ROOT / 'nextflow.config') + ',' + str(project / 'probe.config'),
                       'run', str(project / 'resource_probes.nf'), '-profile', 'local_safe',
                       '-params-file', str(project / 'params.json'), '-work-dir', str(directory / 'work'),
                       '-ansi-log', 'false', '-with-trace', str(directory / 'trace.tsv')]
            run_command(command, directory, 'run')
            with (directory / 'trace.tsv').open(newline='') as handle:
                trace = list(csv.DictReader(handle, delimiter='\t'))
            if len(trace) != len(names) or any(row['status'] != 'COMPLETED' or row['exit'] != '0' for row in trace):
                raise AssertionError(f'{case}: probe trace did not fully complete')
            observations = validate_observations(private, names, case)
            report['probes'][case] = dict(status='PASS', count=len(observations), observations=observations)
            print(f'PASS probes: {case}, {len(observations)} local printf tasks', flush=True)
        report['status'] = 'PASS'
    except Exception as error:
        report['error'] = f'{type(error).__name__}: {error}'
        print(report['error'], file=sys.stderr, flush=True)
    finally:
        report['elapsed_seconds'] = round(time.monotonic() - started, 2)
        write(run / 'report.json', json.dumps(report, indent=2) + '\n')
    return 0 if report['status'] == 'PASS' else 1


if __name__ == '__main__':
    raise SystemExit(main())
