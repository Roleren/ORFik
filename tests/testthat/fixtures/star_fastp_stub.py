#!/usr/bin/env python3
"""Record actual tool arguments and create pipeline outputs; no alignment simulation."""
import json
import os
from pathlib import Path
import shutil
import sys

args = sys.argv[1:]
tool = Path(sys.argv[0]).name
with open(os.environ['ORFIK_TEST_CALLS'], 'a') as stream:
    stream.write(json.dumps({'tool': tool, 'args': args}) + '\n')
if os.environ.get('ORFIK_TEST_FAIL') == tool:
    sys.exit(23)

def values(option):
    if option not in args:
        return []
    result = []
    for value in args[args.index(option) + 1:]:
        if value.startswith('--'):
            break
        result.append(value)
    return result

def one(option):
    return values(option)[0]

if tool == 'fastp':
    for mate in ('1', '2'):
        if '--in' + mate in args:
            shutil.copyfile(one('--in' + mate), one('--out' + mate))
else:
    inputs = values('--readFilesIn')
    for name in inputs:
        assert Path(name).is_file(), name
    prefix = one('--outFileNamePrefix')
    if one('--outReadsUnmapped') == 'Fastx':
        for mate, name in enumerate(inputs, 1):
            shutil.copyfile(name, prefix + 'Unmapped.out.mate' + str(mate))
    if 'BAM' in values('--outSAMtype'):
        suffix = 'Aligned.sortedByCoord.out.bam' if 'SortedByCoordinate' in args else 'Aligned.out.bam'
        Path(prefix + suffix).touch()
    Path(prefix + 'Log.out').write_text('stub log\n')
    Path(prefix + 'SJ.out.tab').touch()
