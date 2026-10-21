#!/usr/bin/env python3
"""make_bundle.py SRC_DIR OUT_DIR : deterministic code bundle for a cluster session that cannot read the shared folder (sent as one pasted message).
Writes OUT_DIR/code.tar.xz, OUT_DIR/msg_kt_code.txt (paste text: header + steps + one base64 data line), OUT_DIR/bundle.sha256.
The tar holds every file of SRC_DIR (not msg_*, not MANIFEST.sha256) plus MANIFEST.sha256 (sha256sum list of those files)."""
import sys, os, io, tarfile, lzma, hashlib, base64
src, out = sys.argv[1:3]; os.makedirs(out, exist_ok=True)
names = sorted(f for f in os.listdir(src) if os.path.isfile(f'{src}/{f}') and not f.startswith('msg_') and f != 'MANIFEST.sha256')
man = ''.join(f"{hashlib.sha256(open(f'{src}/{f}','rb').read()).hexdigest()}  ./{f}\n" for f in names)
open(f'{src}/MANIFEST.sha256', 'w').write(man)
buf = io.BytesIO()
with tarfile.open(fileobj=buf, mode='w', format=tarfile.GNU_FORMAT) as tf:
    for f in names + ['MANIFEST.sha256']:
        data = open(f'{src}/{f}', 'rb').read(); ti = tarfile.TarInfo(f'code/{f}'); ti.size = len(data); ti.mtime = 0; ti.mode = 0o755 if f.endswith(('.sh', '.py', '.sbatch')) else 0o644; ti.uname = ti.gname = ''; ti.uid = ti.gid = 0
        tf.addfile(ti, io.BytesIO(data))
xz = lzma.compress(buf.getvalue(), preset=9 | lzma.PRESET_EXTREME)
open(f'{out}/code.tar.xz', 'wb').write(xz)
b64 = base64.b64encode(xz).decode(); h_tar = hashlib.sha256(xz).hexdigest(); h_line = hashlib.sha256(b64.encode()).hexdigest()
open(f'{out}/bundle.sha256', 'w').write(f'{h_tar}  code.tar.xz\n{h_line}  data line of msg_kt_code.txt\n')
print(len(names), 'files,', len(xz), 'bytes xz,', len(b64), 'base64 chars; tar.xz sha256', h_tar, '; data line sha256', h_line)
open(f'{out}/_data_line.txt', 'w').write(b64)
open(f'{out}/_hashes.txt', 'w').write(f'{h_tar} {h_line}')
