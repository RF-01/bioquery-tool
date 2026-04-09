"""
pipeline.py
End-to-end orchestration: gene symbol -> integrated, summarized result.

Run from the command line for a quick test without Streamlit:
    python pipeline.py TP53
"""

import json
import sys

from fetchers import (
    fetch_uniprot_by_gene,
    parse_uniprot_entry,
    fetch_kegg_pathways,
)
from interpreter import summarize


def run_pipeline(gene_symbol: str, organism_id: int = 9606) -> dict:
    raw = fetch_uniprot_by_gene(gene_symbol, organism_id)
    if raw is None:
        return {"error": f"No reviewed UniProt entry found for gene '{gene_symbol}'"}

    parsed = parse_uniprot_entry(raw)

    # KEGG: UniProt cross-refs already give us 'hsa:7157'-style IDs
    pathways = []
    for kid in parsed["kegg_ids"]:
        try:
            pathways.extend(fetch_kegg_pathways(kid))
        except Exception as e:
            print(f"[warn] KEGG fetch failed for {kid}: {e}", file=sys.stderr)

    # Deduplicate while preserving order
    seen, unique = set(), []
    for p in pathways:
        if p["id"] not in seen:
            seen.add(p["id"])
            unique.append(p)

    return {
        "gene": gene_symbol,
        "parsed": parsed,
        "pathways": unique,
        "summary": summarize(parsed, unique),
    }


if __name__ == "__main__":
    gene = sys.argv[1] if len(sys.argv) > 1 else "TP53"
    result = run_pipeline(gene)

    if "error" in result:
        print(result["error"])
        sys.exit(1)

    s = result["summary"]
    p = result["parsed"]
    print("=" * 70)
    print(f"  {p['protein_name']}  ({', '.join(p['gene_names'])})")
    print(f"  UniProt: {p['accession']}   Length: {p['length']} aa")
    print(f"  GO terms: {s['go_count']}   KEGG pathways: {s['pathway_count']}")
    print("=" * 70)
    print()
    print("INTERPRETIVE SUMMARY")
    print("-" * 70)
    print(s["narrative"])
    print()
    print("TOP GO TERMS BY ASPECT")
    print("-" * 70)
    for aspect, terms in s["ranked_go"].items():
        print(f"\n{aspect}:")
        for t in terms:
            print(f"  {t['id']:<12} {t['term']:<55} [{t['evidence']}]")
    print()
    print("KEGG PATHWAYS")
    print("-" * 70)
    for pw in result["pathways"]:
        print(f"  {pw['id']:<12} {pw['name']}")
