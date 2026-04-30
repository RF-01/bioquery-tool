"""
app.py
Streamlit UI for the unified multi-database query interface.

Run with:
    streamlit run app.py

Phase C additions:
- Handles 'ambiguous' results: shows a selection widget so the user can
  pick the correct protein when a symbol resolves to multiple entries.
- Handles 'not_found' results with a helpful error message.
- Handles 'network' errors with a distinct error message.
- Displays a warning banner for soft issues (gene name mismatch, partial KEGG).
- Uses st.session_state to preserve query results across widget interactions.
"""

import streamlit as st
import altair as alt
import pandas as pd
from pipeline import run_pipeline
from interpreter import clean_function_text, EVIDENCE_TIERS

st.set_page_config(page_title="Unified Gene/Protein Query", layout="wide")

st.title("🧬 Unified Multi-Database Query Interface")
st.caption("Retrieve and interpret functional data from UniProt, GO, and KEGG.")

# ── Session state initialisation ────────────────────────────────────────────
if "result" not in st.session_state:
    st.session_state.result = None
if "last_gene" not in st.session_state:
    st.session_state.last_gene = ""
if "last_organism" not in st.session_state:
    st.session_state.last_organism = 9606

# ── Sidebar ──────────────────────────────────────────────────────────────────
with st.sidebar:
    st.header("Query")
    gene     = st.text_input("Gene symbol", value="TP53")
    organism = st.number_input("Organism (NCBI Taxon ID)", value=9606, step=1)
    go_btn   = st.button("Run query", type="primary")
    st.markdown("---")
    st.markdown(
        "**Examples to try:**\n"
        "- TP53 (tumor suppressor)\n"
        "- BRCA1 (DNA repair)\n"
        "- EGFR (signaling)\n"
        "- INS (insulin)\n"
        "- ALB (serum albumin)\n"
        "- MTOR (non-human: mouse taxon 10090)"
    )

# ── Evidence tier helpers ────────────────────────────────────────────────────
def assign_tier_label(code: str) -> str:
    score = EVIDENCE_TIERS.get((code or "").upper(), 0)
    return {3: "Experimental", 2: "Curated", 1: "Phylogenetic/Sequence", 0: "Electronic (IEA)"}.get(score, "Unknown")

TIER_ORDER  = ["Experimental", "Curated", "Phylogenetic/Sequence", "Electronic (IEA)"]
TIER_COLORS = {
    "Experimental":           "#2ecc71",
    "Curated":                "#3498db",
    "Phylogenetic/Sequence":  "#f39c12",
    "Electronic (IEA)":       "#e74c3c",
}

# ── KEGG category helpers ────────────────────────────────────────────────────
def kegg_category(pathway_id: str) -> str:
    pid = pathway_id.lower()
    if pid.startswith("hsa00") or (pid.startswith("hsa01") and not pid.startswith("hsa015")):
        return "Metabolism"
    elif pid.startswith("hsa015"):
        return "Drug Resistance"
    elif pid.startswith("hsa03"):
        return "Genetic Information Processing"
    elif any(pid.startswith(p) for p in ["hsa040", "hsa041", "hsa042", "hsa043"]):
        return "Cellular Processes & Signaling"
    elif any(pid.startswith(p) for p in ["hsa044", "hsa045", "hsa046", "hsa047", "hsa048", "hsa049"]):
        return "Organismal Systems"
    elif pid.startswith("hsa05"):
        return "Human Diseases"
    else:
        return "Other"

CATEGORY_ORDER  = [
    "Metabolism", "Genetic Information Processing",
    "Cellular Processes & Signaling", "Organismal Systems",
    "Human Diseases", "Drug Resistance", "Other",
]
CATEGORY_COLORS = {
    "Metabolism":                      "#1abc9c",
    "Genetic Information Processing":  "#9b59b6",
    "Cellular Processes & Signaling":  "#3498db",
    "Organismal Systems":              "#f39c12",
    "Human Diseases":                  "#e74c3c",
    "Drug Resistance":                 "#e67e22",
    "Other":                           "#95a5a6",
}

# ── Run query ────────────────────────────────────────────────────────────────
if go_btn:
    st.session_state.last_gene     = gene
    st.session_state.last_organism = int(organism)
    with st.spinner(f"Fetching data for {gene}..."):
        st.session_state.result = run_pipeline(gene, int(organism))

result = st.session_state.result

if result is None:
    st.stop()

# ── Error handling ───────────────────────────────────────────────────────────
error_type = result.get("error_type")

if error_type == "network":
    st.error(f"🌐 **Network error:** {result['error']}")
    st.stop()

if error_type == "not_found":
    st.error(f"🔍 **Gene not found:** {result['error']}")
    st.info(
        "**Tips:**\n"
        "- Check the gene symbol spelling (e.g. `TP53`, not `tp53` or `Tp53`)\n"
        "- Make sure the NCBI Taxon ID matches the organism (human = 9606, mouse = 10090)\n"
        "- Some symbols are aliases — try the official HGNC symbol"
    )
    st.stop()

if error_type == "ambiguous":
    st.warning(
        f"⚠️ **Ambiguous symbol:** `{st.session_state.last_gene}` matches "
        f"{len(result['candidates'])} reviewed UniProt entries. "
        "Please select the protein you meant:"
    )
    options = {
        f"{c['accession']} — {c['protein_name']} ({', '.join(c['gene_names'])})": c["accession"]
        for c in result["candidates"]
    }
    chosen_label = st.selectbox("Select the correct entry:", list(options.keys()))
    if st.button("Use this entry"):
        with st.spinner("Fetching selected entry..."):
            st.session_state.result = run_pipeline(
                st.session_state.last_gene,
                st.session_state.last_organism,
                chosen_accession=options[chosen_label],
            )
        st.rerun()
    st.stop()

# ── Soft warning banner ──────────────────────────────────────────────────────
if "warning" in result:
    st.warning(f"⚠️ {result['warning']}")

# ── Normal result display ────────────────────────────────────────────────────
parsed   = result["parsed"]
summary  = result["summary"]
pathways = result["pathways"]

title_gene = parsed["gene_names"][0] if parsed["gene_names"] else ""
st.subheader(f"{parsed['protein_name']}  ({title_gene})")

c1, c2, c3, c4 = st.columns(4)
c1.metric("UniProt", parsed["accession"])
c2.metric("Length (aa)", parsed["length"])
c3.metric("GO terms", summary["go_count"])
c4.metric("KEGG pathways", summary["pathway_count"])

st.markdown("### 📝 Interpretive Summary")
st.info(summary["narrative"])

tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs([
    "Function", "GO Annotations", "Pathways", "Cross-references",
    "📊 GO Evidence Chart", "🗺️ Pathway Overview",
])

with tab1:
    st.markdown("**Function**")
    for f in parsed["function_text"]:
        st.write("•", clean_function_text(f))
    st.markdown("**Subcellular Location**")
    seen_locs = set()
    for s in parsed["subcellular_location"]:
        if s not in seen_locs:
            seen_locs.add(s)
            st.write("•", s)
    st.markdown("**Keywords**")
    st.write(", ".join(parsed["keywords"][:30]))

with tab2:
    for aspect, terms in summary["ranked_go"].items():
        st.markdown(f"**{aspect}** — top {len(terms)} by evidence quality")
        for t in terms:
            st.write(f"- `{t['id']}` {t['term']}  _(evidence: {t['evidence']})_")

with tab3:
    if pathways:
        for p in pathways:
            st.write(f"- **{p['id']}** — {p['name']}")
    else:
        st.write("No KEGG pathways found.")

with tab4:
    st.write("**KEGG IDs:**", parsed["kegg_ids"])
    st.write("**PDB IDs:**", parsed["pdb_ids"][:20])

# ── Tab 5: GO Evidence Chart ─────────────────────────────────────────────────
with tab5:
    st.markdown("#### GO Term Counts by Evidence Tier")
    st.caption(
        "Shows how many GO annotations fall into each reliability tier. "
        "The ranking rule retains only the top 5 per aspect from the highest tiers, "
        "filtering out the bulk of low-confidence electronic predictions."
    )
    go_anns = parsed.get("go_annotations", [])
    if not go_anns:
        st.write("No GO annotations available.")
    else:
        tier_counts: dict = {t: 0 for t in TIER_ORDER}
        for ann in go_anns:
            label = assign_tier_label(ann.get("evidence", ""))
            tier_counts[label] = tier_counts.get(label, 0) + 1

        df_tier = pd.DataFrame([
            {"Tier": tier, "Count": count}
            for tier, count in tier_counts.items() if count > 0
        ])
        df_tier["Tier"] = pd.Categorical(df_tier["Tier"], categories=TIER_ORDER, ordered=True)
        df_tier = df_tier.sort_values("Tier")
        total = df_tier["Count"].sum()
        df_tier["Percentage"] = (df_tier["Count"] / total * 100).round(1)
        df_tier["Label"] = df_tier.apply(lambda r: f"{r['Count']} ({r['Percentage']}%)", axis=1)

        chart = (
            alt.Chart(df_tier)
            .mark_bar(cornerRadiusTopLeft=4, cornerRadiusTopRight=4)
            .encode(
                x=alt.X("Tier:N", sort=TIER_ORDER, axis=alt.Axis(labelAngle=-20)),
                y=alt.Y("Count:Q", title="Number of GO annotations"),
                color=alt.Color(
                    "Tier:N",
                    scale=alt.Scale(domain=TIER_ORDER, range=[TIER_COLORS[t] for t in TIER_ORDER]),
                    legend=None,
                ),
                tooltip=["Tier", "Count", "Percentage"],
            )
            .properties(height=320)
        )
        text = chart.mark_text(dy=-8, fontSize=13, fontWeight="bold").encode(text="Label:N")
        st.altair_chart(chart + text, use_container_width=True)

        exp_count  = tier_counts.get("Experimental", 0)
        elec_count = tier_counts.get("Electronic (IEA)", 0)
        if total > 0:
            st.info(
                f"**Filtering effect:** Of {total} total GO annotations, "
                f"{exp_count} ({exp_count/total*100:.1f}%) are experimentally supported. "
                f"The ranking rule discards {elec_count} electronic-only predictions "
                f"({elec_count/total*100:.1f}%) when they would otherwise crowd out "
                f"higher-confidence terms."
            )

# ── Tab 6: Pathway Overview ──────────────────────────────────────────────────
with tab6:
    st.markdown("#### KEGG Pathway Category Overview")
    st.caption(
        "Groups all retrieved KEGG pathways by their top-level biological class. "
        "Provides a quick view of which biological domains this protein participates in."
    )
    if not pathways:
        st.write("No KEGG pathways found for this gene.")
    else:
        cat_counts:   dict = {}
        cat_pathways: dict = {}
        for p in pathways:
            cat = kegg_category(p["id"])
            cat_counts[cat]   = cat_counts.get(cat, 0) + 1
            cat_pathways.setdefault(cat, []).append(p)

        df_cat = pd.DataFrame([
            {"Category": cat, "Count": count}
            for cat, count in cat_counts.items()
        ])
        df_cat["Category"] = pd.Categorical(
            df_cat["Category"], categories=CATEGORY_ORDER, ordered=True
        )
        df_cat = df_cat.sort_values("Category")

        chart2 = (
            alt.Chart(df_cat)
            .mark_bar(cornerRadiusTopLeft=4, cornerRadiusTopRight=4)
            .encode(
                x=alt.X("Category:N", sort=CATEGORY_ORDER, axis=alt.Axis(labelAngle=-25)),
                y=alt.Y("Count:Q", title="Number of KEGG pathways"),
                color=alt.Color(
                    "Category:N",
                    scale=alt.Scale(
                        domain=CATEGORY_ORDER,
                        range=[CATEGORY_COLORS[c] for c in CATEGORY_ORDER],
                    ),
                    legend=None,
                ),
                tooltip=["Category", "Count"],
            )
            .properties(height=320)
        )
        text2 = chart2.mark_text(dy=-8, fontSize=13, fontWeight="bold").encode(text="Count:Q")
        st.altair_chart(chart2 + text2, use_container_width=True)

        st.markdown("---")
        st.markdown("**Pathways by category** _(click to expand)_")
        for cat in CATEGORY_ORDER:
            if cat in cat_pathways:
                with st.expander(f"{cat} — {cat_counts[cat]} pathway(s)"):
                    for p in cat_pathways[cat]:
                        name = p["name"].split(" - ")[0]
                        st.write(f"- **{p['id']}** — {name}")
