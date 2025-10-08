import io
import os
import re
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import requests
from Bio.PDB import PDBIO, PDBParser

# Defaults
ALPHAFOLD_PDB_URL_TEMPLATE = "https://alphafold.ebi.ac.uk/files/{chain_id}.pdb"
DEFAULT_CACHE_DIR_ENV = "CATH_VIS_CACHE_DIR"
DEFAULT_CACHE_DIR = Path.home() / ".cache" / "cath-vis"


def get_cache_dir() -> Path:
    env_path = os.getenv(DEFAULT_CACHE_DIR_ENV)
    if env_path:
        return Path(env_path).expanduser().resolve()
    return DEFAULT_CACHE_DIR


def ensure_cache_dir(cache_dir: Path | None = None) -> Path:
    cache_path = cache_dir or get_cache_dir()
    cache_path.mkdir(parents=True, exist_ok=True)
    return cache_path


def normalize_uniprot_to_chain_id(uniprot_or_chain_id: str) -> str:
    # Accepts either a full AF chain_id or a bare UniProt accession; returns AF chain_id
    if uniprot_or_chain_id.startswith("AF-") and uniprot_or_chain_id.endswith("-model_v4"):
        return uniprot_or_chain_id
    acc = uniprot_or_chain_id.strip()
    if not acc:
        raise ValueError("Empty UniProt accession/chain_id provided")
    # AlphaFold canonical monomer file naming pattern
    return f"AF-{acc}-F1-model_v4"


def parse_boundaries(boundary_string: str) -> List[Tuple[int, int]]:
    if not isinstance(boundary_string, str) or boundary_string.strip() == "" or boundary_string.lower() == "nan":
        return []
    ranges: List[Tuple[int, int]] = []
    # segments separated by '_' and within each segment commas separate pairs like '5-42,60-120'
    for segment in boundary_string.split("_"):
        for pair in segment.split(","):
            pair = pair.strip()
            if not pair:
                continue
            m = re.match(r"^(\d+)-(\d+)$", pair)
            if not m:
                continue
            start, end = int(m.group(1)), int(m.group(2))
            if end < start:
                start, end = end, start
            ranges.append((start, end))
    # merge overlapping/adjacent ranges for cleanliness
    if not ranges:
        return []
    ranges.sort()
    merged: List[Tuple[int, int]] = []
    cur_start, cur_end = ranges[0]
    for s, e in ranges[1:]:
        if s <= cur_end + 1:
            cur_end = max(cur_end, e)
        else:
            merged.append((cur_start, cur_end))
            cur_start, cur_end = s, e
    merged.append((cur_start, cur_end))
    return merged


def residue_color_map(
    structure,
    hmm_boundaries: Sequence[Tuple[int, int]] | None,
    foldseek_boundaries: Sequence[Tuple[int, int]] | None,
    colors: Dict[str, str] | None = None,
) -> Dict[int, str]:
    color_palette = {
        "overlap": "#776bcd",  # soft purple
        "hmm": "#22a7f0",      # soft blue
        "ted": "#e14b31",      # soft red
        "none": "#E0E0E0",     # light gray
    }
    if colors:
        color_palette.update(colors)

    hmm = list(hmm_boundaries or [])
    ted = list(foldseek_boundaries or [])

    def in_ranges(resnum: int, ranges: Iterable[Tuple[int, int]]) -> bool:
        for start, end in ranges:
            if start <= resnum <= end:
                return True
        return False

    mapped: Dict[int, str] = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] == " ":  # standard amino acid
                    resnum = residue.id[1]
                    is_hmm = in_ranges(resnum, hmm)
                    is_ted = in_ranges(resnum, ted)
                    if is_hmm and is_ted:
                        mapped[resnum] = color_palette["overlap"]
                    elif is_hmm:
                        mapped[resnum] = color_palette["hmm"]
                    elif is_ted:
                        mapped[resnum] = color_palette["ted"]
                    else:
                        mapped[resnum] = color_palette["none"]
    return mapped


def fetch_alphafold_pdb(chain_id: str, cache_dir: Path | None = None) -> str:
    cache_dir = ensure_cache_dir(cache_dir)
    cache_file = cache_dir / f"{chain_id}.pdb"
    if cache_file.exists() and cache_file.stat().st_size > 0:
        return cache_file.read_text()

    url = ALPHAFOLD_PDB_URL_TEMPLATE.format(chain_id=chain_id)
    resp = requests.get(url, timeout=60)
    if resp.status_code != 200:
        raise ValueError(f"Failed to download structure for {chain_id}: HTTP {resp.status_code}")
    cache_file.write_text(resp.text)
    return resp.text


def load_structure_from_pdb_string(pdb_text: str):
    parser = PDBParser(QUIET=True)
    return parser.get_structure("protein", io.StringIO(pdb_text))


def structure_to_pdb_string(structure) -> str:
    pdb_io = io.StringIO()
    writer = PDBIO()
    writer.set_structure(structure)
    writer.save(pdb_io)
    return pdb_io.getvalue()


def build_py3dmol_viewer(structure, color_by_residue: Dict[int, str], width: int = 800, height: int = 600):
    import py3Dmol  # imported lazily to avoid heavy import cost when unused

    viewer = py3Dmol.view(width=width, height=height)
    pdb_string = structure_to_pdb_string(structure)
    viewer.addModel(pdb_string, "pdb")
    viewer.setStyle({"cartoon": {"colorscheme": {"prop": "resi", "map": color_by_residue}}})
    viewer.zoomTo()
    return viewer
