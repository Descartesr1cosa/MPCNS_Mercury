#!/usr/bin/env python3
"""Run fresh/restart CT initialization on two tiny blocks (one or two MPI ranks)."""
import argparse
from pathlib import Path
import re
import shutil
import struct
import math
import subprocess
import tempfile


def prepare(root, ranks, rotated=False):
    app = Path(__file__).resolve().parents[1]
    setup = root / 'CASE/setup'
    shutil.copytree(app / '9999setup', setup)
    general = setup / '0_General_Setting.txt'
    text = general.read_text()
    for key, value in [('Read_Method', 0), ('max_Nstep', 0), ('corner_treat', 0)]:
        text = re.sub(rf'({key}\s+)\d+', rf'\g<1>{value}', text)
    general.write_text(text)
    control = setup / '3_MHD_Control.txt'
    control.write_text(re.sub(r'(is_Crustal_Bfield\s+)0', r'\g<1>1', control.read_text()))
    g = root / 'CASE/geometry'
    bc = g / 'boundary_condition'
    bc.mkdir(parents=True)
    (bc / 'boundary_name.fvbnd').write_text('FVBND 1 3\nNo_Boundary_Condition\nFarfield\nBOUNDARIES\n')
    for rank in range(ranks):
        blocks = list(range(2)) if ranks == 1 else [rank]
        with (g / f'grid{rank:4d}.grd').open('w') as f:
            f.write(f'{len(blocks)}\n' + 'Fluid\n' * len(blocks) + '6 6 6\n' * len(blocks))
            for b in blocks:
                def point(i, j, k):
                    if b == 1 and rotated:
                        return (2.1+0.02*j, 0.3-0.02*i, 0.2+0.02*k)
                    return (2+0.1*b+0.02*i, 0.2+0.02*j, 0.2+0.02*k)
                for d in range(3):
                    f.write(' '.join(str(point(i, j, k)[d])
                                     for k in range(6) for j in range(6) for i in range(6)) + '\n')
        with (bc / f'boundary{rank:4d}.txt').open('w') as f:
            for b in blocks:
                f.write('5\n')
                for d in range(3):
                    for side in [0, 5]:
                        normal = 1 if b == 1 and rotated else 0
                        if d == normal and side == (5 if b == 0 else 0):
                            continue
                        bounds = [0, 5, 0, 5, 0, 5]
                        bounds[2*d:2*d+2] = [side, side]
                        f.write(' '.join(map(str, bounds)) + ' 2\n')
        if ranks == 1:
            receiver = '6 1 1 1 -1 -6' if rotated else '1 1 1 6 -1 -6'
            (bc / f'inner{rank:4d}.txt').write_text(
                f'2\n1 1\nheader\n6 6 1 6 -1 -6 1 Fluid\n{receiver} 2 Fluid\n'
                f'{receiver} 2 Fluid\n6 6 1 6 -1 -6 1 Fluid\n')
            (bc / f'parallel{rank:4d}.txt').write_text('2 0\n0 0\nheader\n')
        else:
            (bc / f'inner{rank:4d}.txt').write_text('0\n0\nheader\n')
            side = 6 if rank == 0 else 1
            face = '6 1 1 1 -1 -6' if rank == 1 and rotated else f'{side} {side} 1 6 -1 -6'
            (bc / f'parallel{rank:4d}.txt').write_text(
                f'1 1\n1\nheader\n{face} {1-rank} {2*(rank+1)} {2*(2-rank)} Fluid\n')


def header_offset(data):
    if data[:8] != b'MPCNSRST':
        raise RuntimeError('bad restart magic')
    version, = struct.unpack_from('=i', data, 8)
    if version != 2:
        raise RuntimeError('bad restart version')
    length, = struct.unpack_from('=i', data, 32)
    if data[36:36+length] != b'Mars_Bind_Badd_IMF_Dipoles_v1':
        raise RuntimeError('wrong magnetic convention')
    return 36+length


def records(data):
    blocks, count = struct.unpack_from('=2i', data, 24)
    offset = header_offset(data)
    for _ in range(count):
        start = offset
        length, = struct.unpack_from('=i', data, offset)
        offset += 4
        name = data[offset:offset+length].decode()
        offset += length
        _, components, ghost = struct.unpack_from('=3i', data, offset)
        offset += 12
        slots = []
        for block in range(blocks):
            li, lj, lk, hi, hj, hk, active = struct.unpack_from('=7i', data, offset)
            offset += 28
            if not active:
                continue
            for i in range(li, hi):
                for j in range(lj, hj):
                    for k in range(lk, hk):
                        for c in range(components):
                            inside = 0 <= i < hi-ghost and 0 <= j < hj-ghost and 0 <= k < hk-ghost
                            slots.append((offset, inside))
                            offset += 8
        yield name, start, offset, slots
    if offset != len(data):
        raise RuntimeError('unexpected trailing bytes')


def fields(data):
    return {name: [struct.unpack_from('=d', data, off)[0] for off, inner in slots if inner]
            for name, _, _, slots in records(data)}


def magnetic_ranges(output):
    line = next(s for s in output.splitlines() if 'B    x=[' in s)
    return [float(v) for v in re.findall(r'[-+]?\d+\.\d+e[-+]\d+', line)]


def check_replacement(before, after):
    # Setup changes only Badd. Stored/evolved Bind and ion state stay unchanged.
    for filename, data in before.items():
        old, new = fields(data), fields(after[filename])
        for name in old:
            if (name.startswith('U_') or name in ['B_xi', 'B_eta', 'B_zeta']) and old[name] != new[name]:
                raise RuntimeError('background replacement changed induced field/ion state')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    parser.add_argument('--ranks', type=int, choices=[1, 2], default=2)
    parser.add_argument('--rotated', action='store_true', help='rotate the second block axes by 90 degrees')
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='mars-dipoles-smoke-') as tmp:
        root = Path(tmp)
        prepare(root, args.ranks, args.rotated)
        command = ['mpiexec', '-n', str(args.ranks), str(args.executable.resolve())]
        def run(expected_cache=None, expected_failure=None):
            try:
                result = subprocess.run(command, cwd=root, capture_output=True, text=True, timeout=120)
            except subprocess.TimeoutExpired as error:
                raise RuntimeError(f'MPI test timed out: {error.stdout!r}\n{error.stderr!r}') from error
            if expected_failure:
                if result.returncode == 0 or expected_failure not in result.stdout+result.stderr:
                    raise RuntimeError('expected rejection missing: '+result.stdout+result.stderr)
                return
            if result.returncode:
                raise RuntimeError(result.stdout + result.stderr)
            if expected_cache and '[Mars Dipole cache] '+expected_cache not in result.stdout:
                raise RuntimeError('unexpected cache path: '+result.stdout)
            line = next(s for s in result.stdout.splitlines() if '[Mars Badd]' in s)
            print(line)
            run.output = result.stdout
            return line
        if 'fresh Bind=0' not in run():
            raise RuntimeError('fresh Badd initialization missing')
        match = re.search(r'\|J\|max=([0-9.eE+-]+)', run.output)
        if not match or float(match[1]) != 0:
            raise RuntimeError('prescribed background produced a solved current')
        files = list((root / 'DATA').glob('flow_field*.bin'))
        if len(files) != args.ranks:
            raise RuntimeError('expected one restart field file per rank')
        before = {f.name: f.read_bytes() for f in files}
        for data in before.values():
            if set(fields(data)) != {'U_'+name for name in ['Ar','CO2','CO','C','H','N2','N','NO','O','O2']} | {'B_xi','B_eta','B_zeta','Badd_xi','Badd_eta','Badd_zeta'}:
                raise RuntimeError('checkpoint must contain ions, induced flux and background flux')
            for name in ['B_xi', 'B_eta', 'B_zeta']:
                if any(v != 0 for v in fields(data)[name]):
                    raise RuntimeError('fresh induced field contains background')
        general = root / 'CASE/setup/0_General_Setting.txt'
        general.write_text(re.sub(r'(continue_calc\s+)0', r'\g<1>1', general.read_text()))
        control = root / 'CASE/setup/3_MHD_Control.txt'
        control.write_text(re.sub(r'(Read_Badd_Dipoles\s+)0', r'\g<1>1', control.read_text()))
        if 'restart Bind retained' not in run('hit'):
            raise RuntimeError('restart Badd initialization missing')
        after = {f.name: f.read_bytes() for f in (root / 'DATA').glob('flow_field*.bin')}
        if before != after:
            raise RuntimeError('restart changed conserved state/face flux (possible double addition)')
        print(f'{args.ranks} MPI rank(s), two blocks, 4840 sources: restart fields byte-identical')
        def snapshot():
            return {f.name: f.read_bytes() for f in (root / 'DATA').glob('flow_field*.bin')}
        initial_ranges = magnetic_ranges(run.output)
        # Replace the file with a physically different (twice as strong) model.
        dipoles = root / 'CASE/setup/Dipoles'
        rows = dipoles.read_text().splitlines()
        for n, row in enumerate(rows):
            if row.startswith('#') or row.startswith('MARS_') or not row.strip():
                continue
            values = list(map(float, row.split()))
            values[3:6] = [2*v for v in values[3:6]]
            rows[n] = ' '.join(format(v, '.16e') for v in values)
        dipoles.write_text('\n'.join(rows)+'\n')
        run('miss')
        changed = snapshot()
        check_replacement(after, changed)
        if magnetic_ranges(run.output) == initial_ranges:
            raise RuntimeError('changed model had no effect')
        control = root / 'CASE/setup/3_MHD_Control.txt'
        control.write_text(re.sub(r'(is_Crustal_Bfield\s+)1', r'\g<1>0', control.read_text()))
        run()
        disabled = snapshot()
        check_replacement(changed, disabled)
        # Total cell field must equal the latest IMF, once, in Cartesian axes.
        disabled_ranges = magnetic_ranges(run.output)
        if disabled_ranges != [-0.56,-0.56,0.83,0.83,0.0,0.0]:
            raise RuntimeError('IMF missing/doubled when dipoles disabled')
        control.write_text(re.sub(r'(is_Crustal_Bfield\s+)0', r'\g<1>1', control.read_text()))
        run('hit')
        check_replacement(disabled, snapshot())
        print('Changed model, disable on restart, and re-enable on restart: passed')
        enabled = snapshot()
        # Stored background is a postprocessing snapshot, never a restart input.
        for name,data in enabled.items():
            damaged = bytearray(data)
            for field,_,_,slots in records(data):
                if field.startswith('Badd_'):
                    for offset,_ in slots:
                        struct.pack_into('=d',damaged,offset,1e90)
            (root / 'DATA' / name).write_bytes(damaged)
        run('hit')
        if snapshot() != enabled:
            raise RuntimeError('restart used the saved background instead of current setup')
        # This is a breaking release: reject old untagged files and incomplete
        # tagged schemas instead of silently interpreting old total B as Bind.
        for mode,error in [('version','unsupported version'),('missing','missing required field')]:
            for name,data in enabled.items():
                if mode == 'version':
                    modified = data[:8]+struct.pack('=i',1)+data[12:32]+data[header_offset(data):]
                else:
                    kept = [data[start:end] for field,start,end,_ in records(data) if field != 'Badd_xi']
                    modified = data[:28]+struct.pack('=i',len(kept))+data[32:header_offset(data)]+b''.join(kept)
                (root / 'DATA' / name).write_bytes(modified)
            run(expected_failure=error)
        for name,data in enabled.items():
            (root / 'DATA' / name).write_bytes(data)
        print('Saved Badd ignored; obsolete/incomplete checkpoint formats rejected: passed')
        stable = snapshot()
        cache = root / 'DATA/Badd_Dipoles/dipole0000.bin'
        damaged = bytearray(cache.read_bytes())
        damaged[-1] ^= 1
        cache.write_bytes(damaged)
        run('miss')
        if stable != snapshot():
            raise RuntimeError('recomputing corrupt cache changed same-model checkpoint')
        cache.unlink()
        run('miss')
        if stable != snapshot():
            raise RuntimeError('recomputing missing cache changed checkpoint')
        run('hit')
        print('Cache hit, changed model invalidation, corrupt and missing rank cache fallback: passed')
        stable_ranges = magnetic_ranges(run.output)
        # IMF is rebuilt independently of the dipole cache, exactly once.
        initial = root / 'CASE/setup/1_Initial_condition.txt'
        original_initial = initial.read_text()
        initial.write_text(original_initial.replace('-1.68E-9', '-3.36E-9').replace('2.49E-9', '4.98E-9'))
        run('hit')
        new_imf = snapshot()
        check_replacement(stable, new_imf)
        delta = [-0.56,-0.56,0.83,0.83,0.0,0.0]
        if any(not math.isclose(y-x,z,abs_tol=0.002) for x,y,z in
               zip(stable_ranges,magnetic_ranges(run.output),delta)):
            raise RuntimeError('IMF change not applied exactly once')
        for filename in stable:
            old,new,imf = fields(stable[filename]),fields(new_imf[filename]),fields(disabled[filename])
            for axis in ['xi','eta','zeta']:
                name = 'Badd_'+axis
                if any(not math.isclose(y-x,z,abs_tol=1e-14) for x,y,z in zip(old[name],new[name],imf[name])):
                    raise RuntimeError('saved Badd does not reflect current IMF exactly once')
        initial.write_text(original_initial)
        # Exercise the physical EMF: moving plasma across nonuniform prescribed
        # dipoles must generate Bind even though its initial value/current is zero.
        settings = re.sub(r'(continue_calc\s+)1', r'\g<1>0', general.read_text())
        settings = re.sub(r'(max_Nstep\s+)0', r'\g<1>1', settings)
        settings = re.sub(r'(CFL\s+)0.5', r'\g<1>0.0001', settings)
        settings = re.sub(r'(output_step\s+)1000', r'\g<1>1', settings)
        general.write_text(settings)
        run('hit')
        evolved = snapshot()
        if any(struct.unpack_from('=i', data, 12)[0] != 1 for data in evolved.values()):
            raise RuntimeError('one-step checkpoint was not written')
        if max(abs(v) for data in evolved.values() for name in ['B_xi','B_eta','B_zeta'] for v in fields(data)[name]) < 1e-14:
            raise RuntimeError('ideal EMF failed to use the prescribed background')
        general.write_text(re.sub(r'(continue_calc\s+)0', r'\g<1>1', general.read_text()))
        initial.write_text(original_initial.replace('-1.68E-9', '-3.36E-9'))
        control.write_text(re.sub(r'(is_Crustal_Bfield\s+)1', r'\g<1>0', control.read_text()))
        run()
        check_replacement(evolved, snapshot())
        print('IMF-only cache reuse, one-step induction, and nonzero Bind retained on background change: passed')
        # Uniform stationary plasma and zero Bind: Hall/background jumps must
        # not diffuse the prescribed nonuniform dipoles into the evolved field.
        general.write_text(re.sub(r'(continue_calc\s+)1', r'\g<1>0', general.read_text()))
        control.write_text(re.sub(r'(is_Crustal_Bfield\s+)0', r'\g<1>1', control.read_text()))
        initial.write_text(original_initial.replace('5.0E5', '0.0'))
        boundary = root / 'CASE/setup/2_Boundary_Condition.txt'
        boundary.write_text(boundary.read_text().replace('5.0E5', '0.0'))
        run('hit')
        resting = snapshot()
        if any(not math.isfinite(v) or abs(v) > 1e-14 for data in resting.values()
               for name in ['B_xi','B_eta','B_zeta'] for v in fields(data)[name]):
            raise RuntimeError('stationary background generated spurious induction')
        print('Stationary nonuniform background with Hall enabled: no spurious induction')




if __name__ == '__main__':
    main()
