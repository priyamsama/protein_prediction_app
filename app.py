# app.py
import streamlit as st
import requests
import py3Dmol

st.title("Protein 3D Structure Prediction with ESMFold")

sequence = st.text_area("Enter protein sequence in FASTA format (without header):")

if st.button("Predict Structure"):
    if not sequence.strip():
        st.error("Please enter a valid protein sequence.")
    else:
        st.info("Sending sequence to ESMFold API...")
        response = requests.post(
            "https://api.esmatlas.com/foldSequence/v1/pdb/",
            headers={"Content-Type": "application/x-www-form-urlencoded"},
            data=sequence,
        )
        if response.status_code == 200:
            pdb_data = response.text
            viewer = py3Dmol.view(width=600, height=400)
            viewer.addModel(pdb_data, "pdb")
            viewer.setStyle({'cartoon': {'color': 'spectrum'}})
            viewer.zoomTo()
            viewer.show()
            st.components.v1.html(viewer._make_html(), height=400)
        else:
            st.error("Prediction failed. Please try again.")
