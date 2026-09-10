#!/usr/bin/env python3
"""Run fresh/restart CT initialization on two tiny blocks (one or two MPI ranks)."""
import argparse
from pathlib import Path
import re
import shutil
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
        def run():
            try:
                result = subprocess.run(command, cwd=root, capture_output=True, text=True, timeout=120)
            except subprocess.TimeoutExpired as error:
                raise RuntimeError(f'MPI test timed out: {error.stdout!r}\n{error.stderr!r}') from error
            if result.returncode:
                raise RuntimeError(result.stdout + result.stderr)
            line = next(s for s in result.stdout.splitlines() if '[Mars Badd]' in s)
            print(line)
            return line
        if 'added to initial total field' not in run():
            raise RuntimeError('fresh Badd initialization missing')
        files = list((root / 'DATA').glob('flow_field*.bin'))
        if len(files) != args.ranks:
            raise RuntimeError('expected one restart field file per rank')
        before = {f.name: f.read_bytes() for f in files}
        general = root / 'CASE/setup/0_General_Setting.txt'
        general.write_text(re.sub(r'(continue_calc\s+)0', r'\g<1>1', general.read_text()))
        if 'restart total field retained' not in run():
            raise RuntimeError('restart Badd initialization missing')
        after = {f.name: f.read_bytes() for f in (root / 'DATA').glob('flow_field*.bin')}
        if before != after:
            raise RuntimeError('restart changed conserved state/face flux (possible double addition)')
        print(f'{args.ranks} MPI rank(s), two blocks, 4840 sources: restart fields byte-identical')


if __name__ == '__main__':
    main()
