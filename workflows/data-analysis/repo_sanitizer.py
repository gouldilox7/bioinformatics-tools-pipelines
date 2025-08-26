#!/usr/bin/env python3
"""
repo_sanitizer.py — scan & sanitize a repository for names, paths, IDs, and secrets,
plus normalize *file references* to 'study' + extension (.csv, .tsv, .txt, .job, .sh, .batch, .R, .py).

Official docs referenced:
- argparse: https://docs.python.org/3/library/argparse.html
- re: https://docs.python.org/3/library/re.html
- json: https://docs.python.org/3/library/json.html
- PyYAML safe_load: https://pyyaml.org/wiki/PyYAMLDocumentation
- Jupyter nbformat (schema reference for clearing outputs): https://nbformat.readthedocs.io/en/latest/format_description.html
"""
import argparse, re, sys, os, json
from typing import List, Dict, Any, Tuple

try:
    import yaml  # PyYAML
except Exception:
    yaml = None

DEFAULT_EXCLUDE_DIRS = {
    '.git','.hg','.svn','__pycache__','node_modules',
    '.mypy_cache','.ruff_cache','.pytest_cache','.venv','venv','env',
    '.Rproj.user','.ipynb_checkpoints'
}

DEFAULT_TEXT_EXT = {
    '.py','.r','.R','.Rmd','.md','.txt','.csv','.tsv','.json','.yml','.yaml',
    '.toml','.ini','.cfg','.sh','.bash','.zsh','.ps1','.bat','.sql','.sv','.xml',
    '.html','.htm','.css','.js','.ts','.ipynb','.pl','.jl','.m','.c','.cpp','.h','.hpp'
}

SECRET_PATTERNS = [
    ("AWS Access Key", r"AKIA[0-9A-Z]{16}"),
    ("AWS Secret", r"(?i)aws(.{0,20})?(secret|access).{0,20}?['\"][A-Za-z0-9/+=]{40}['\"]"),
    ("Generic API Key", r"(?i)(api[_-]?key|token|secret|password)\\s*[:=]\\s*['\"][A-Za-z0-9_\\-]{16,}['\"]"),
    ("Slack Token", r"xox[abpr]-[A-Za-z0-9-]{10,48}"),
    ("GitHub Token", r"gh[pous]_[A-Za-z0-9]{36,}"),
    ("Google API Key", r"AIza[0-9A-Za-z\\-_]{35}"),
    ("Private Key Block", r"-----BEGIN (?:RSA|EC|OPENSSH|PGP) PRIVATE KEY-----"),
]

def is_probably_text(path: str, max_check=4096) -> bool:
    try:
        with open(path, 'rb') as f:
            chunk = f.read(max_check)
        if not chunk:
            return True
        if b'\\x00' in chunk:
            return False
        chunk.decode('utf-8', errors='ignore')
        return True
    except Exception:
        return False

def load_config(cfg_path: str) -> Dict[str, Any]:
    if yaml is None:
        raise RuntimeError("PyYAML is not installed. Install with: pip install pyyaml")
    with open(cfg_path, 'r', encoding='utf-8') as f:
        cfg = yaml.safe_load(f) or {}
    # Defaults
    cfg.setdefault('exclude_dirs', [])
    cfg.setdefault('include_ext', [])
    cfg.setdefault('exclude_ext', [])
    cfg.setdefault('replacements', [])
    cfg.setdefault('size_limit_mb', 5)
    cfg.setdefault('scan_only_patterns', [])
    # New: normalize file references
    nfr = cfg.get('normalize_file_refs') or {}
    nfr.setdefault('enabled', False)
    nfr.setdefault('extensions', [])
    nfr.setdefault('rename_files', False)  # content-only by default
    cfg['normalize_file_refs'] = nfr
    return cfg

def compile_rules(replacements):
    compiled = []
    for i, r in enumerate(replacements):
        pat = r.get('pattern')
        repl = r.get('replacement', '<REDACTED>')
        flags = 0
        for f in (r.get('flags') or []):
            lf = f.lower()
            if lf == 'i': flags |= re.IGNORECASE
            elif lf == 'm': flags |= re.MULTILINE
            elif lf == 's': flags |= re.DOTALL
        try:
            compiled.append((re.compile(pat, flags), repl, r.get('label', f'rule_{i}')))
        except re.error as e:
            print(f"[WARN] Invalid regex in replacements: {pat} — {e}", file=sys.stderr)
    return compiled

def compile_scan_only(scan_patterns):
    compiled = []
    for i, p in enumerate(scan_patterns):
        try:
            flags = 0
            for f in (p.get('flags') or []):
                lf = f.lower()
                if lf == 'i': flags |= re.IGNORECASE
                elif lf == 'm': flags |= re.MULTILINE
                elif lf == 's': flags |= re.DOTALL
            compiled.append({'label': p.get('label', f'scan_{i}'), 're': re.compile(p['pattern'], flags)})
        except re.error as e:
            print(f"[WARN] Invalid scan-only regex: {p.get('pattern')} — {e}", file=sys.stderr)
    return compiled

def should_skip_file(path: str, cfg) -> bool:
    ext = os.path.splitext(path)[1]
    size_limit = int(cfg.get('size_limit_mb', 5)) * 1024 * 1024
    try:
        if os.path.getsize(path) > size_limit:
            return True
    except Exception:
        return True
    if ext and cfg.get('exclude_ext'):
        if ext.lower() in {e.lower() for e in cfg['exclude_ext']}:
            return True
    if cfg.get('include_ext'):
        if ext.lower() not in {e.lower() for e in cfg['include_ext']}:
            return True
    return not is_probably_text(path)

def iter_files(root: str, cfg):
    exclude_dirs = DEFAULT_EXCLUDE_DIRS.union(set(cfg.get('exclude_dirs', [])))
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if d not in exclude_dirs]
        for fn in filenames:
            yield os.path.join(dirpath, fn)

def build_file_ref_regex(exts: List[str]) -> re.Pattern:
    # exts like ['.csv','.tsv','.R','.py'] -> 'csv|tsv|R|py' (escape dots)
    exts_clean = [re.escape(e.lstrip('.')) for e in exts]
    ext_alt = '|'.join(exts_clean)
    # Match file references with optional directories. We preserve dir path and replace only the last component's stem.
    # Boundaries: start or whitespace or quote; end or whitespace/quote.
    pattern = rf'(?P<prefix>(?<=^)|(?<=[\s\"\'\(\)=:]))' \
              rf'(?P<dir>(?:[A-Za-z]:)?(?:[^\s\"\'<>|]*[\\/])*)' \
              rf'(?P<name>[^\\/\s\"\'<>|]+?)' \
              rf'(?P<ext>\.(?:{ext_alt}))' \
              rf'(?P<suffix>(?=$)|(?=[\s\"\'<>|,)\]]))'
    return re.compile(pattern)

def normalize_file_refs(data: str, exts: List[str]) -> Tuple[str, int, List[Tuple[str,str]]]:
    rx = build_file_ref_regex(exts)
    changes = []
    def _sub(m: re.Match):
        before = m.group(0)
        new = f"{m.group('prefix')}{m.group('dir')}study{m.group('ext')}{m.group('suffix')}"
        if before != new:
            changes.append((before, new))
        return new
    new_data, n = rx.subn(_sub, data)
    return new_data, n, changes

def scan_data_string(data: str, compiled_rules, scan_only_res, nfr_cfg):
    findings = []
    # Built-in secret detection
    for name, pat in SECRET_PATTERNS:
        for m in re.finditer(pat, data):
            snippet = data[max(0, m.start()-20):m.end()+20].replace('\\n',' ')
            findings.append({'type': 'secret', 'name': name, 'span': (m.start(), m.end()), 'snippet': snippet})
    # Report-only patterns
    for s in scan_only_res:
        for m in s['re'].finditer(data):
            snippet = data[max(0, m.start()-20):m.end()+20].replace('\\n',' ')
            findings.append({'type': 'pattern', 'name': s['label'], 'span': (m.start(), m.end()), 'snippet': snippet})
    # Replacements (report)
    for rx, repl, label in compiled_rules:
        for m in rx.finditer(data):
            s = data[m.start():m.end()]
            findings.append({'type': 'replace', 'name': label, 'span': (m.start(), m.end()), 'before': s[:100], 'after': rx.sub(repl, s)[:100]})
    # Normalization report (no change yet)
    if nfr_cfg.get('enabled') and nfr_cfg.get('extensions'):
        _rx = build_file_ref_regex(nfr_cfg['extensions'])
        for m in _rx.finditer(data):
            seg = m.group(0)
            # Skip if it already matches 'study.<ext>'
            if m.group('name') != 'study':
                proposal = f"{m.group('dir')}study{m.group('ext')}"
                findings.append({'type': 'normalize_ref', 'name': 'file_ref_to_study', 'span': (m.start(), m.end()), 'before': seg[:120], 'after': proposal[:120]})
    return findings

def apply_replacements_and_normalize(data: str, compiled_rules, nfr_cfg) -> Tuple[str, int, int, List[Tuple[str,str]]]:
    total_repl = 0
    total_norm = 0
    # Apply configured replacements
    for rx, repl, _ in compiled_rules:
        data, n = rx.subn(repl, data)
        total_repl += n
    # Normalize file references to 'study.<ext>'
    changes = []
    if nfr_cfg.get('enabled') and nfr_cfg.get('extensions'):
        data, n, changes = normalize_file_refs(data, nfr_cfg['extensions'])
        total_norm += n
    return data, total_repl, total_norm, changes

def strip_notebook_outputs(json_str: str) -> str:
    try:
        nb = json.loads(json_str)
    except Exception:
        return json_str
    changed = False
    if isinstance(nb, dict) and isinstance(nb.get('cells'), list):
        for cell in nb['cells']:
            if isinstance(cell, dict) and cell.get('cell_type') == 'code':
                if 'outputs' in cell and cell['outputs']:
                    cell['outputs'] = []
                    changed = True
                if 'execution_count' in cell and cell['execution_count'] is not None:
                    cell['execution_count'] = None
                    changed = True
    if changed:
        return json.dumps(nb, ensure_ascii=False, indent=1)
    return json_str

def process_file(path: str, compiled_rules, scan_only_res, cfg, apply: bool, strip_ipynb: bool):
    ext = os.path.splitext(path)[1].lower()
    findings = []
    changed = 0
    try:
        with open(path, 'r', encoding='utf-8', errors='ignore') as f:
            data = f.read()
    except Exception:
        return findings, changed

    findings.extend(scan_data_string(data, compiled_rules, scan_only_res, cfg['normalize_file_refs']))

    if apply:
        new_data, n_repl, n_norm, changes = apply_replacements_and_normalize(data, compiled_rules, cfg['normalize_file_refs'])
        if ext == '.ipynb' and strip_ipynb:
            new_data = strip_notebook_outputs(new_data)
        if new_data != data:
            with open(path, 'w', encoding='utf-8') as f:
                f.write(new_data)
            changed = 1
            if n_norm:
                # add a summarized finding for normalization actions
                findings.append({'type': 'normalized_count', 'name': 'file_ref_to_study', 'count': n_norm})
    return findings, changed

def main():
    ap = argparse.ArgumentParser(description="Scan and sanitize a repository; normalize file references to 'study.<ext>'.")
    ap.add_argument('--root', required=True, help='Root directory to scan')
    ap.add_argument('--config', required=True, help='YAML config file (patterns.yml)')
    mode = ap.add_mutually_exclusive_group(required=False)
    mode.add_argument('--dry-run', action='store_true', help='Scan only (default)')
    mode.add_argument('--apply', action='store_true', help='Apply replacements/normalization in-place')
    ap.add_argument('--strip-notebooks', action='store_true', help='Strip outputs/execution_count from .ipynb')
    ap.add_argument('--report', default='sanitize_report.json', help='Path for JSON report')
    args = ap.parse_args()

    root = args.root
    if not os.path.isdir(root):
        print(f"[ERR] Root not found: {root}", file=sys.stderr)
        sys.exit(2)

    if yaml is None:
        print("[ERR] PyYAML is required. Install with: python -m pip install pyyaml", file=sys.stderr)
        sys.exit(3)

    cfg = load_config(args.config)
    compiled_rules = compile_rules(cfg.get('replacements', []))
    scan_only_res = compile_scan_only(cfg.get('scan_only_patterns', []))

    scanned = 0
    changed_files = 0
    all_findings = []

    for path in iter_files(root, cfg):
        if should_skip_file(path, cfg):
            continue
        scanned += 1
        findings, changed = process_file(path, compiled_rules, scan_only_res, cfg, args.apply, args.strip_notebooks)
        if findings:
            all_findings.append({'file': os.path.abspath(path), 'findings': findings})
        changed_files += changed

    with open(args.report, 'w', encoding='utf-8') as f:
        json.dump({'root': os.path.abspath(root),'scanned_files': scanned,'changed_files': changed_files,'findings': all_findings}, f, indent=2)

    mode_text = "DRY-RUN" if (args.dry_run or not args.apply) else "APPLY"
    print(f"[{mode_text}] Scanned {scanned} files; modified {changed_files}.")
    print(f"[REPORT] {args.report}")
    if cfg.get('normalize_file_refs',{}).get('enabled'):
        print("File-ref normalization is ENABLED. (Content-only; no on-disk renames.)")

if __name__ == '__main__':
    main()
