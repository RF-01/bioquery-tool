"""
pipeline.py
End-to-end orchestration: gene symbol -> integrated, summarized result.

Phase C additions:
- Detects ambiguous gene symbols (multiple reviewed UniProt hits) and returns
  a structured 'ambiguous' result so the UI can ask the user to pick one.
- Detects not-found / obsolete symbols and returns a clear error.
- Wraps KEGG calls in per-gene error handling so a KEGG failure never
  crashes the whole pipeline.
- Accepts an optional accession argument so the UI can re-run with a
  specific UniProt entry after resolving ambiguity.

Run from the command line for a quick test without Streamlit:
    python pipeline.py TP53
    python pipeline.py ALB          # will warn about ambiguity, uses first hit
"""

import json
import sys

from fetchers import (
    fetch_uniprot_by_gene,
    parse_uniprot_entry,
    fetch_kegg_pathways,
    FetchError,
)
from interpreter import summarize


def run_pipeline(
    gene_symbol: str,
    organism_id: int = 9606,
    chosen_accession: Optional[str] = None,
) -> dict:
    """
    Run the full pipeline for *gene_symbol* in *organism_id*.

    Return dict keys:
      On success:   gene, parsed, pathways, summary
                    + optional 'warning' (str) for soft issues
      On ambiguity: error_type='ambiguous', candidates=[{accession, protein_name, gene_names}]
      On not found: error_type='not_found',  error=<message>
      On network:   error_type='network',    error=<message>
    """
    # ── 1. Fetch reviewed UniProt hits ──────────────────────────────────────
    try:
        hits = fetch_uniprot_by_gene(gene_symbol, organism_id, max_hits=5)
    except FetchError as exc:
        return {"error": str(exc), "error_type": "network"}

    # Not found / obsolete symbol
    if not hits:
        return {
            "error": (
                f"No reviewed UniProt entry found for gene symbol '{gene_symbol}' "
                f"in organism {organism_id}. "
                "The symbol may be obsolete, misspelled, or absent from Swiss-Prot."
            ),
            "error_type": "not_found",
        }

    # Ambiguity: multiple hits and the caller hasn't chosen one yet
    if len(hits) > 1 and chosen_accession is None:
        candidates = []
        for h in hits:
            acc  = h.get("primaryAccession", "?")
            prot = (
                h.get("proteinDescription", {})
                 .get("recommendedName", {})
                 .get("fullName", {})
                 .get("value", "Unknown protein")
            )
            genes = [
                g.get("geneName", {}).get("value", "")
                for g in h.get("genes", [])
                if g.get("geneName", {}).get("value")
            ]
            candidates.append({"accession": acc, "protein_name": prot, "gene_names": genes})
        return {"error_type": "ambiguous", "candidates": candidates}

    # Pick the entry: either the caller's chosen accession or the top hit
    if chosen_accession:
        entry = next((h for h in hits if h.get("primaryAccession") == chosen_accession), hits[0])
    else:
        entry = hits[0]

    parsed = parse_uniprot_entry(entry)

    # Warn if the primary gene name doesn't match what was queried
    # (catches cases like querying "ALB" and getting a different gene back)
    primary_gene = (parsed["gene_names"] or [""])[0].upper()
    warning = None
    if primary_gene and primary_gene != gene_symbol.upper():
        warning = (
            f"The top UniProt hit for '{gene_symbol}' is annotated as gene "
            f"'{primary_gene}' ({parsed['accession']}). "
            "This may indicate an ambiguous or alias symbol. "
            "Consider specifying the UniProt accession directly if this is unexpected."
        )

    # ── 2. KEGG pathway lookup ───────────────────────────────────────────────
    pathways      = []
    kegg_warnings = []
    for kid in parsed["kegg_ids"]:
        try:
            pathways.extend(fetch_kegg_pathways(kid))
        except FetchError as exc:
            kegg_warnings.append(f"KEGG fetch failed for {kid}: {exc}")
        except Exception as exc:
            kegg_warnings.append(f"Unexpected error for {kid}: {exc}")

    # Deduplicate while preserving order
    seen, unique = set(), []
    for p in pathways:
        if p["id"] not in seen:
            seen.add(p["id"])
            unique.append(p)

    # Merge any KEGG warnings into the main warning string
    if kegg_warnings:
        kegg_note = "KEGG data may be incomplete: " + "; ".join(kegg_warnings)
        warning = (warning + " | " + kegg_note) if warning else kegg_note

    result = {
        "gene":     gene_symbol,
        "parsed":   parsed,
        "pathways": unique,
        "summary":  summarize(parsed, unique),
    }
    if warning:
        result["warning"] = warning
    return result


# Allow `chosen_accession` type hint without importing Optional at top
from typing import Optional
run_pipeline.__annotations__["chosen_accession"] = Optional[str]


if __name__ == "__main__":
    gene = sys.argv[1] if len(sys.argv) > 1 else "TP53"
    result = run_pipeline(gene)

    if result.get("error_type") == "ambiguous":
        print(f"Ambiguous symbol '{gene}' — multiple reviewed entries found:")
        for c in result["candidates"]:
            print(f"  {c['accession']}  {c['protein_name']}  ({', '.join(c['gene_names'])})")
        print("Using the first hit. Pass a specific accession to avoid this.")
        # Re-run with the first candidate
        result = run_pipeline(gene, chosen_accession=result["candidates"][0]["accession"])

    if "error" in result:
        print(result["error"])
        sys.exit(1)

    if "warning" in result:
        print(f"[warning] {result['warning']}\n")

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
