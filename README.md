# Unified Multi-Database Query Interface

A bioinformatics tool that retrieves gene/protein information from **UniProt**,
**Gene Ontology (GO)**, and **KEGG** in a single query, then produces a
rule-based interpretive summary instead of just dumping raw records.

## Setup

```bash
pip install -r requirements.txt
```

## Quick test (no UI)

```bash
python pipeline.py TP53
```

This prints a formatted summary to the terminal — useful for verifying the
APIs work and for capturing screenshots for the progress report.

Try other genes too:

```bash
python pipeline.py BRCA1
python pipeline.py EGFR
```

## Streamlit UI

```bash
streamlit run app.py
```

Then open the URL it prints (usually http://localhost:8501).

## Files

| File             | Purpose                                                    |
|------------------|------------------------------------------------------------|
| `fetchers.py`    | UniProt search/parse, KEGG pathway lookup                  |
| `interpreter.py` | Rule-based GO ranking + narrative composition              |
| `pipeline.py`    | Orchestrator + CLI entry point                             |
| `app.py`         | Streamlit UI                                               |

## How the interpretation works

GO annotations are grouped by aspect (Biological Process, Molecular Function,
Cellular Component) and ranked using evidence-code tiers:

- **Tier 3** — experimental (IDA, IMP, IPI, EXP, …)
- **Tier 2** — author statement / curator (TAS, NAS, IC)
- **Tier 1** — phylogenetic / sequence-based (ISS, IBA, …)
- **Tier 0** — electronic only (IEA)

This filters out noisy electronic annotations and surfaces the
biology that humans have actually verified — the core "added value"
beyond a plain database fetch.
