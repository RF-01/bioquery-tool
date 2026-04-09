"""
fetchers.py
API clients for UniProt and KEGG.
GO annotations are extracted from UniProt's cross-references, which contain
the GO ID, term name, aspect (P/F/C), and evidence code — avoiding a second API.
"""

import requests
from typing import Optional, Dict, Any, List

UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb"
KEGG_BASE = "https://rest.kegg.jp"
HEADERS = {"Accept": "application/json"}
TIMEOUT = 30


# ---------- UniProt ----------

def fetch_uniprot_by_gene(gene_symbol: str, organism_id: int = 9606) -> Optional[Dict[str, Any]]:
    """
    Search UniProt for a reviewed (Swiss-Prot) entry matching gene symbol + organism.
    Default organism is human (taxon 9606). Returns the top hit as a dict, or None.
    """
    query = f"(gene:{gene_symbol}) AND (organism_id:{organism_id}) AND (reviewed:true)"
    params = {"query": query, "format": "json", "size": 1}
    r = requests.get(f"{UNIPROT_BASE}/search", params=params, headers=HEADERS, timeout=TIMEOUT)
    r.raise_for_status()
    results = r.json().get("results", [])
    return results[0] if results else None


def parse_uniprot_entry(entry: Dict[str, Any]) -> Dict[str, Any]:
    """
    Flatten a raw UniProt JSON entry into a clean dict the rest of the pipeline uses.
    """
    out: Dict[str, Any] = {
        "accession": entry.get("primaryAccession"),
        "uniprot_id": entry.get("uniProtkbId"),
        "protein_name": None,
        "gene_names": [],
        "organism": None,
        "length": entry.get("sequence", {}).get("length"),
        "function_text": [],
        "subcellular_location": [],
        "keywords": [],
        "go_annotations": [],   # list of {id, term, aspect, evidence}
        "kegg_ids": [],
        "pdb_ids": [],
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
        db = ref.get("database")
        rid = ref.get("id")
        if db == "GO":
            term, aspect, evidence = None, None, None
            for p in ref.get("properties", []):
                k, v = p.get("key"), p.get("value")
                if k == "GoTerm" and v:
                    # Format "P:DNA repair" — first char is aspect code
                    if ":" in v:
                        code, term_name = v.split(":", 1)
                        aspect = {"P": "Biological Process",
                                  "F": "Molecular Function",
                                  "C": "Cellular Component"}.get(code, code)
                        term = term_name
                    else:
                        term = v
                elif k == "GoEvidenceType" and v:
                    evidence = v.split(":", 1)[0]   # e.g. "IDA:UniProtKB" -> "IDA"
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
    """
    # Step 1: get pathway IDs linked to this gene
    r = requests.get(f"{KEGG_BASE}/link/pathway/{kegg_gene_id}", timeout=TIMEOUT)
    r.raise_for_status()
    pathway_ids = []
    for line in r.text.strip().splitlines():
        parts = line.split("\t")
        if len(parts) == 2:
            pid = parts[1].replace("path:", "")
            pathway_ids.append(pid)
    if not pathway_ids:
        return []

    # Step 2: pull the master human pathway list once and look up names
    org = kegg_gene_id.split(":")[0]   # e.g. 'hsa'
    r2 = requests.get(f"{KEGG_BASE}/list/pathway/{org}", timeout=TIMEOUT)
    r2.raise_for_status()
    name_map = {}
    for line in r2.text.strip().splitlines():
        parts = line.split("\t")
        if len(parts) == 2:
            key = parts[0].replace("path:", "")
            name_map[key] = parts[1]

    return [{"id": pid, "name": name_map.get(pid, pid)} for pid in pathway_ids]
