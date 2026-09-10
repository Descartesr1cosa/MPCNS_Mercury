#!/usr/bin/env python3
"""Convert the supplied Langlais M23/-20/14, DIP_FORM=1 solution to Mars SI dipoles."""
import argparse
import hashlib
import math
from pathlib import Path

RADIUS_KM = 3393.5  # R_M in the author's forward.c; independent of solver L_ref.


def convert(source, destination):
    raw = source.read_bytes()
    lines = raw.decode('utf-8-sig').splitlines()
    start = None
    for i, line in enumerate(lines):
        words = line.split()
        if len(words) != 3:
            continue
        try:
            spacing, count, depth = float(words[0]), int(words[1]), float(words[2])
        except ValueError:
            continue
        if count == 4840 and spacing == 2.9189 and depth == 20.0:
            start = i
            break
    if start is None:
        raise ValueError('Expected M23/-20/14 header: 2.9189 4840 20.000')
    rows = [list(map(float, line.split())) for line in lines[start+1:] if line.strip()]
    if len(rows) != count:
        raise ValueError(f'Expected {count} records, found {len(rows)}')
    result = [
        '# Langlais, Purucker & Mandea (2004), DOI:10.1029/2003JE002048, M23/-20/14.',
        f'# Source: {source.name}; SHA256={hashlib.sha256(raw).hexdigest()}',
        '# Converted with tools/convert_langlais2004.py using the supplied forward.c, DIP_FORM=1.',
        '# Mars body-fixed Cartesian axes: X=(lat 0, east lon 0), Y=(lat 0, east lon 90), Z=north.',
        '# No body-to-solar/simulation rotation applied. Rotate positions AND moments if required.',
        '# r=(3393.5-20) km; V=(3393.5*pi/180*2.9189)^2*thickness_km [km^3].',
        '# m=V*1e9*(Mr*e_r+Mt*e_theta+Mp*e_phi) [A m^2]; e_theta=south, e_phi=east.',
        '# x_m y_m z_m mx_Am2 my_Am2 mz_Am2 core_m; mu0/(4*pi) NOT included.',
        f'MARS_DIPOLES_V1 {count}',
    ]
    volumes = []
    for number, row in enumerate(rows, start+2):
        if len(row) != 6 or not all(math.isfinite(x) for x in row):
            raise ValueError(f'{source}:{number}: expected six finite values')
        lon, lat, thickness, mr, mt, mp = row
        if not (0 <= lon < 360 and -90 <= lat <= 90 and thickness > 0):
            raise ValueError(f'{source}:{number}: invalid coordinates/thickness')
        theta, phi = math.radians(90-lat), math.radians(lon)
        st, ct, sp, cp = math.sin(theta), math.cos(theta), math.sin(phi), math.cos(phi)
        er, et, ep = (st*cp, st*sp, ct), (ct*cp, ct*sp, -st), (-sp, cp, 0.0)
        volume = (RADIUS_KM*math.pi/180*spacing)**2 * thickness
        volumes.append(volume)
        position = [(RADIUS_KM-depth)*1000*x for x in er]
        moment = [volume*1e9*(mr*er[d]+mt*et[d]+mp*ep[d]) for d in range(3)]
        result.append(' '.join(f'{x:.16e}' for x in position+moment+[0.0]))
    if source.resolve() == destination.resolve():
        raise ValueError('Keep the source separate from the converted Dipoles')
    destination.write_text('\n'.join(result)+'\n')
    print(f'{count} dipoles -> {destination}; source radius={(RADIUS_KM-depth):.1f} km; '
          f'volume range={min(volumes):.12g}..{max(volumes):.12g} km^3')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('destination', type=Path)
    args = parser.parse_args()
    convert(args.source, args.destination)
