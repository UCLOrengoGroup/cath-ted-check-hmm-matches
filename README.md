# CATH/TED domain visualization (lean)

Minimal tool to visualize HMM vs TED domain boundaries on AlphaFold structures. Analyses live under `analyses/`.

## Setup

```bash
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

Optional cache dir for AlphaFold PDB files (defaults to ~/.cache/cath-vis):

```bash
export CATH_VIS_CACHE_DIR=/path/to/cache
```

## Usage in notebooks

```python
import pandas as pd
from domain_vis import visualize_domains_and_plot

df = pd.read_csv("results/domain_comparison_results_human.processed.annotated.tsv", sep="\t")
visualize_domains_and_plot(df)
```

- Color toggle: pLDDT or packing density; if a `match` column exists, a Match/No Match mode is available.
- Search accepts UniProt accessions or full AlphaFold IDs; accessions normalize to `AF-<ACC>-F1-model_v4`.

## Modules

- `domain_vis/core.py`: cache-aware AlphaFold fetch, boundary parsing/merging, viewer builder
- `domain_vis/visualize.py`: unified interactive viewer (exported via `domain_vis/__init__.py`)

## Analyses (optional)

Scripts are kept separate in `analyses/` to keep the viewer lean:

```bash
python analyses/plot_density.py --input data/domain_comparison_results_human.tsv --output density_human.png
python analyses/scatterplot.py --input ted-evalue.txt --outdir results/plots --ymax 80
python analyses/migrate_upgrade_v2.py
python analyses/migrate_annotate.py
```

## Data inputs

- `data/domain_comparison_results_human.tsv`
- `results/domain_comparison_results_human.processed.tsv`
- `results/domain_comparison_results_human.processed.annotated.tsv`

## Notes

- AlphaFold PDBs are cached; delete cache files to refresh.
- Boundary strings like `5-42,60-120_200-240` are parsed and merged.



