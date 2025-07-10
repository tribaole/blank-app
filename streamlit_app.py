import streamlit as st
import requests
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction

# --- Utility functions ---
def get_spliceai_prediction(seq, pos, ref, alt):
    spliceai_url = "https://spliceailookup-api.broadinstitute.org/spliceai"
    payload = {
        "sequence": seq,
        "position": pos,
        "ref": ref,
        "alt": alt
    }
    try:
        response = requests.post(spliceai_url, json=payload)
        if response.status_code == 200:
            return response.json()
        else:
            return {"error": f"SpliceAI request failed: {response.status_code}"}
    except Exception as e:
        return {"error": str(e)}


def generate_candidate_oligos(seq, length=20):
    return [seq[i:i + length] for i in range(len(seq) - length + 1)]


def calculate_gc(seq):
    return round(gc_fraction(seq) * 100, 2)

# --- Streamlit App ---
st.title("Splice-Switching Antisense Oligonucleotide Designer")

st.markdown("""
Upload a target pre-mRNA sequence, and select a region of interest to design candidate SSOs.
This demo integrates **SpliceAI** to help predict functional impact of binding.
""")

# User input
sequence = st.text_area("Paste target pre-mRNA sequence", height=200)

region_start = st.number_input("Start position of target region (1-based)", min_value=1, value=1)
region_end = st.number_input("End position of target region (1-based)", min_value=1, value=30)

if st.button("Generate SSOs"):
    if not sequence or region_start > region_end or region_end > len(sequence):
        st.error("Please provide a valid sequence and region.")
    else:
        region = sequence[region_start - 1:region_end]
        st.write(f"Selected region ({region_start}-{region_end}):")
        st.code(region)

        candidates = generate_candidate_oligos(region)
        st.markdown("### Candidate SSO Sequences")

        for idx, oligo in enumerate(candidates):
            gc = calculate_gc(oligo)
            st.markdown(f"**SSO #{idx + 1}**")
            st.code(oligo)
            st.write(f"GC content: {gc}%")

            # Call SpliceAI (placeholder mutation at middle of oligo)
            middle = len(oligo) // 2
            ref_base = oligo[middle]
            alt_base = 'A' if ref_base != 'A' else 'G'
            spliceai_result = get_spliceai_prediction(oligo, middle, ref_base, alt_base)

            if "error" in spliceai_result:
                st.warning(spliceai_result['error'])
            else:
                st.json(spliceai_result)

        st.success(f"Generated {len(candidates)} candidate SSOs.")
