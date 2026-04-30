"""
fetchers.py
API clients for UniProt and KEGG.
GO annotations are extracted from UniProt's cross-references, which contain
the GO ID, term name, aspect (P/F/C), and evidence code — avoiding a second API.

Phase C additions:
- fetch_uniprot_by_gene now returns ALL reviewed hits (up to 5) so the caller
  can detect ambiguous symbols and let the user choose.
- Network errors are caught and re-raised as FetchError with a human-readable message.
- Obsolete / not-found symbols return an empty list rather than None.
"""

import requests
from typing import Optional, Dict, Any, List

UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb"
KEGG_BASE    = "https://rest.kegg.jp"
HEADERS      = {"Accept": "application/json"}
TIMEOUT      = 30


class FetchError(RuntimeError):
    """Raised when a network call fails after retries."""


# ---------- UniProt ----------

def fetch_uniprot_by_gene(
    gene_symbol: str,
    organism_id: int = 9606,
    max_hits: int = 5,
) -> List[Dict[str, Any]]:
    """
    Search UniProt for reviewed (Swiss-Prot) entries matching gene symbol + organism.
    Returns up to *max_hits* results as a list of raw entry dicts.

    Returning a list (instead of a single entry) lets the pipeline:
      - return an error when nothing is found (empty list)
      - warn the user when several proteins match the same symbol (ambiguous)
      - just use results[0] in the normal single-match case
    """
    query  = f"(gene:{gene_symbol}) AND (organism_id:{organism_id}) AND (reviewed:true)"
    params = {"query": query, "format": "json", "size": max_hits}
    try:
        r = requests.get(
            f"{UNIPROT_BASE}/search", params=params, headers=HEADERS, timeout=TIMEOUT
        )
        r.raise_for_status()
    except requests.exceptions.ConnectionError:
        raise FetchError(
            "Could not reach UniProt. Please check your internet connection."
        )
    except requests.exceptions.Timeout:
        raise FetchError(
            "UniProt request timed out. The server may be busy — try again shortly."
        )
    except requests.exceptions.HTTPError as exc:
        raise FetchError(f"UniProt returned an error: {exc.response.status_code}")

    return r.json().get("results", [])


def parse_uniprot_entry(entry: Dict[str, Any]) -> Dict[str, Any]:
    """
    Flatten a raw UniProt JSON entry into a clean dict the rest of the pipeline uses.
    """
    out: Dict[str, Any] = {
        "accession":           entry.get("primaryAccession"),
        "uniprot_id":          entry.get("uniProtkbId"),
        "protein_name":        None,
        "gene_names":          [],
        "organism":            None,
        "length":              entry.get("sequence", {}).get("length"),
        "function_text":       [],
        "subcellular_location":[],
        "keywords":            [],
        "go_annotations":      [],   # list of {id, term, aspect, evidence}
        "kegg_ids":            [],
        "pdb_ids":             [],
    }

    # Protein name
    rec = entry.get("proteinDescription", {}).get("recommendedName", {})
    if rec:
        out["protein_name"] = rec.get("fullName", {}).get("value")

    # Gene names
    for g in entry.get("genes", []):
        nm = g.get("geneName", {}).get("value")
        if nm:
            out["gene_names"].append(nm)

    # Organism
    out["organism"] = entry.get("organism", {}).get("scientificName")

    # Function + Subcellular location from comments
    for c in entry.get("comments", []):
        ctype = c.get("commentType")
        if ctype == "FUNCTION":
            for t in c.get("texts", []):
                if t.get("value"):
                    out["function_text"].append(t["value"])
        elif ctype == "SUBCELLULAR LOCATION":
            for sl in c.get("subcellularLocations", []):
                loc = sl.get("location", {}).get("value")
                if loc:
                    out["subcellular_location"].append(loc)

    # Keywords
    for k in entry.get("keywords", []):
        if k.get("name"):
            out["keywords"].append(k["name"])

    # Cross references — GO, KEGG, PDB
    for ref in entry.get("uniProtKBCrossReferences", []):
        db  = ref.get("database")
        rid = ref.get("id")
        if db == "GO":
            term, aspect, evidence = None, None, None
            for p in ref.get("properties", []):
                k, v = p.get("key"), p.get("value")
                if k == "GoTerm" and v:
                    if ":" in v:
                        code, term_name = v.split(":", 1)
                        aspect = {"P": "Biological Process",
                                  "F": "Molecular Function",
                                  "C": "Cellular Component"}.get(code, code)
                        term = term_name
                    else:
                        term = v
                elif k == "GoEvidenceType" and v:
                    evidence = v.split(":", 1)[0]
            out["go_annotations"].append(
                {"id": rid, "term": term, "aspect": aspect, "evidence": evidence}
            )
        elif db == "KEGG":
            out["kegg_ids"].append(rid)
        elif db == "PDB":
            out["pdb_ids"].append(rid)

    return out


# ---------- KEGG ----------

def fetch_kegg_pathways(kegg_gene_id: str) -> List[Dict[str, str]]:
    """
    Given a KEGG gene id like 'hsa:7157', return all linked pathways as
    a list of {id, name} dicts.
    Raises FetchError on network problems.
    """
    try:
        r = requests.get(
            f"{KEGG_BASE}/link/pathway/{kegg_gene_id}", timeout=TIMEOUT
        )
        r.raise_for_status()
    except requests.exceptions.ConnectionError:
        raise FetchError(
            "Could not reach KEGG. Please check your internet connection."
        )
    except requests.exceptions.Timeout:
        raise FetchError("KEGG request timed out — try again shortly.")
    except requests.exceptions.HTTPError as exc:
        raise FetchError(f"KEGG returned an error: {exc.response.status_code}")

    pathway_ids = []
    for line in r.text.strip().splitlines():
        parts = line.split("\t")
        if len(parts) == 2:
            pid = parts[1].replace("path:", "")
            pathway_ids.append(pid)
    if not pathway_ids:
        return []

    org = kegg_gene_id.split(":")[0]
    try:
        r2 = requests.get(
            f"{KEGG_BASE}/list/pathway/{org}", timeout=TIMEOUT
        )
        r2.raise_for_status()
    except Exception:
        # Non-fatal: return IDs without names rather than crashing
        return [{"id": pid, "name": pid} for pid in pathway_ids]

    name_map = {}
    for line in r2.text.strip().splitlines():
        parts = line.split("\t")
        if len(parts) == 2:
            key = parts[0].replace("path:", "")
            name_map[key] = parts[1]

    return [{"id": pid, "name": name_map.get(pid, pid)} for pid in pathway_ids]
