"""
interpreter.py
Rule-based interpretation layer.

Strategy:
- Group GO annotations by aspect (Biological Process / Molecular Function /
  Cellular Component).
- Rank within each aspect using GO evidence-code tiers
  (experimental > curated > phylogenetic/sequence > electronic).
- Compose a narrative paragraph from the parsed UniProt data + ranked terms +
  KEGG pathways.

This is intentionally rule-based — no LLM. The "interpretation" is the act of
filtering noise (low-evidence IEA terms), prioritizing experimentally supported
biology, and structuring the result into a readable summary.
"""

from collections import defaultdict
from typing import Dict, Any, List

# GO evidence-code quality tiers (higher = stronger biological support)
EVIDENCE_TIERS = {
    # Experimental
    "EXP": 3, "IDA": 3, "IPI": 3, "IMP": 3, "IGI": 3, "IEP": 3,
    "HTP": 3, "HDA": 3, "HMP": 3, "HGI": 3, "HEP": 3,
    # Author statement / curator
    "TAS": 2, "NAS": 2, "IC":  2,
    # Phylogenetic / sequence / computational analysis
    "ISS": 1, "ISO": 1, "ISA": 1, "ISM": 1, "IGC": 1,
    "IBA": 1, "IBD": 1, "IKR": 1, "IRD": 1, "RCA": 1,
    # Electronic (lowest)
    "IEA": 0,
}


def evidence_score(code: str) -> int:
    return EVIDENCE_TIERS.get((code or "").upper(), 0)


def rank_go_terms(go_annotations: List[Dict[str, Any]],
                  top_n: int = 5) -> Dict[str, List[Dict]]:
    """
    Group GO annotations by aspect and return the top N per aspect,
    ranked by evidence quality (descending), then term name.
    """
    by_aspect = defaultdict(list)
    for ann in go_annotations:
        by_aspect[ann.get("aspect") or "Unknown"].append(ann)

    ranked = {}
    for aspect, anns in by_aspect.items():
        ranked[aspect] = sorted(
            anns,
            key=lambda a: (-evidence_score(a.get("evidence")), a.get("term") or "")
        )[:top_n]
    return ranked


def summarize(parsed: Dict[str, Any],
              pathways: List[Dict[str, str]]) -> Dict[str, Any]:
    """
    Build a structured summary: a narrative paragraph plus ranked GO terms
    and pathway counts. Consumed by the UI and the CLI.
    """
    ranked_go = rank_go_terms(parsed.get("go_annotations", []), top_n=5)

    parts = []
    name = parsed.get("protein_name") or "This protein"
    gene = (parsed.get("gene_names") or [None])[0] or parsed.get("accession")
    organism = parsed.get("organism") or "the organism"
    length = parsed.get("length", "?")
    parts.append(f"{name} ({gene}, {organism}) is a protein of {length} amino acids.")

    if parsed.get("function_text"):
        first = parsed["function_text"][0]
        if len(first) > 350:
            first = first[:347] + "..."
        parts.append(f"Function: {first}")

    if parsed.get("subcellular_location"):
        locs = ", ".join(parsed["subcellular_location"][:4])
        parts.append(f"Localized to: {locs}.")

    bp = ranked_go.get("Biological Process", [])
    if bp:
        parts.append("Top biological processes: " +
                     ", ".join(t["term"] for t in bp) + ".")

    mf = ranked_go.get("Molecular Function", [])
    if mf:
        parts.append("Top molecular functions: " +
                     ", ".join(t["term"] for t in mf) + ".")

    cc = ranked_go.get("Cellular Component", [])
    if cc:
        parts.append("Cellular components: " +
                     ", ".join(t["term"] for t in cc) + ".")

    if pathways:
        # KEGG names often look like "p53 signaling pathway - Homo sapiens (human)"
        names = [p["name"].split(" - ")[0] for p in pathways[:6]]
        parts.append(f"Participates in KEGG pathways including: {', '.join(names)}.")

    return {
        "narrative": " ".join(parts),
        "ranked_go": ranked_go,
        "pathway_count": len(pathways),
        "go_count": len(parsed.get("go_annotations", [])),
    }
