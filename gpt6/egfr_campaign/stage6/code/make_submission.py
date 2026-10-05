#!/usr/bin/env python3
"""Create a <=20-entry CSV from the ranked FASTA, skipping official rejections.
This offline formatter does not run, predict, or bypass Adaptyv's novelty checks.
"""
from __future__ import annotations
import argparse
import csv
from pathlib import Path
import re
import sys

def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name: str | None = None
    parts: list[str] = []
    for number, raw in enumerate(path.read_text(encoding='utf-8').splitlines(), 1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith('>'):
            if name is not None:
                records.append((name, ''.join(parts)))
            name, parts = line[1:], []
            if not re.fullmatch(r'[A-Za-z0-9_]+', name):
                raise ValueError(f'Invalid FASTA identifier on line {number}: {name!r}')
        else:
            if name is None:
                raise ValueError(f'Sequence before first header on line {number}')
            parts.append(line)
    if name is not None:
        records.append((name, ''.join(parts)))
    if not records:
        raise ValueError('The FASTA contains no records.')
    if len({n for n, _ in records}) != len(records):
        raise ValueError('Duplicate identifiers in the FASTA.')
    if len({s for _, s in records}) != len(records):
        raise ValueError('Duplicate sequences in the FASTA.')
    for n, s in records:
        if not 10 <= len(s) <= 250 or not set(s) <= set('ACDEFGHIKLMNPQRSTVWY'):
            raise ValueError(f'Invalid amino acids or length in {n}.')
    return records

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fasta', type=Path, default=Path('egfr_ranked_100.fasta'))
    parser.add_argument('--reject-file', type=Path, help='One rejected public identifier per line; blank lines and # comments are ignored.')
    parser.add_argument('--count', type=int, choices=range(1, 21), default=20, metavar='1..20')
    parser.add_argument('--output', type=Path, default=Path('egfr_track3_revised_top20.csv'))
    args = parser.parse_args()
    try:
        records = read_fasta(args.fasta)
        rejected: set[str] = set()
        if args.reject_file is not None:
            rejected = {s.strip() for s in args.reject_file.read_text(encoding='utf-8').splitlines() if s.strip() and not s.lstrip().startswith('#')}
        unknown = rejected - {n for n, _ in records}
        if unknown:
            raise ValueError('Unknown rejected identifiers: ' + ', '.join(sorted(unknown)))
        chosen = [(n, s) for n, s in records if n not in rejected][:args.count]
        if len(chosen) < args.count:
            raise ValueError(f'Only {len(chosen)} non-rejected entries remain; requested {args.count}.')
        if args.output.resolve() in {args.fasta.resolve(), *( [args.reject_file.resolve()] if args.reject_file else [])}:
            raise ValueError('The output must not overwrite an input.')
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open('w', newline='', encoding='utf-8') as stream:
            writer = csv.writer(stream, lineterminator='\n')
            writer.writerow(['name', 'sequence', 'molecule_class'])
            writer.writerows((n, s, 'nanobody') for n, s in chosen)
        print(f'Wrote {len(chosen)} entries to {args.output}. Official eligibility and biological function remain unverified.')
    except (OSError, ValueError) as error:
        parser.exit(2, f'Error: {error}\n')

if __name__ == '__main__':
    main()
