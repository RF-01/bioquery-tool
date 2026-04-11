"""
app.py
Streamlit UI for the unified multi-database query interface.

Run with:
    streamlit run app.py
"""

import streamlit as st
from pipeline import run_pipeline
from interpreter import clean_function_text

st.set_page_config(page_title="Unified Gene/Protein Query", layout="wide")

st.title("🧬 Unified Multi-Database Query Interface")
st.caption("Retrieve and interpret functional data from UniProt, GO, and KEGG.")

with st.sidebar:
    st.header("Query")
    gene = st.text_input("Gene symbol", value="TP53")
    organism = st.number_input("Organism (NCBI Taxon ID)", value=9606, step=1)
    go_btn = st.button("Run query", type="primary")
    st.markdown("---")
    st.markdown(
        "**Examples to try:**\n"
        "- TP53 (tumor suppressor)\n"
        "- BRCA1 (DNA repair)\n"
        "- EGFR (signaling)\n"
        "- INS (insulin)"
    )

if go_btn:
    with st.spinner(f"Fetching data for {gene}..."):
        result = run_pipeline(gene, int(organism))

    if "error" in result:
        st.error(result["error"])
    else:
        parsed = result["parsed"]
        summary = result["summary"]
        pathways = result["pathways"]

        # Header card
        title_gene = parsed["gene_names"][0] if parsed["gene_names"] else ""
        st.subheader(f"{parsed['protein_name']}  ({title_gene})")

        c1, c2, c3, c4 = st.columns(4)
        c1.metric("UniProt", parsed["accession"])
        c2.metric("Length (aa)", parsed["length"])
        c3.metric("GO terms", summary["go_count"])
        c4.metric("KEGG pathways", summary["pathway_count"])

        # Narrative summary (the "interpretation")
        st.markdown("### 📝 Interpretive Summary")
        st.info(summary["narrative"])

        # Tabs for the underlying integrated data
        tab1, tab2, tab3, tab4 = st.tabs(
            ["Function", "GO Annotations", "Pathways", "Cross-references"]
        )

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
                    st.write(
                        f"- `{t['id']}` {t['term']}  "
                        f"_(evidence: {t['evidence']})_"
                    )

        with tab3:
            if pathways:
                for p in pathways:
                    st.write(f"- **{p['id']}** — {p['name']}")
            else:
                st.write("No KEGG pathways found.")

        with tab4:
            st.write("**KEGG IDs:**", parsed["kegg_ids"])
            st.write("**PDB IDs:**", parsed["pdb_ids"][:20])
