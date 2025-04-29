import streamlit as st
import numpy as np
import py3Dmol
import requests
import time
import io
from stmol import showmol
import biotite.structure.io.pdb as pdb
from PIL import Image

# Set page configuration
st.set_page_config(
    page_title="Protein Structure Prediction App",
    page_icon="🧬",
    layout="wide",
)

# Custom CSS
st.markdown("""
<style>
    .main-header {
        font-size: 2.5rem;
        color: #4257B2;
        margin-bottom: 1rem;
    }
    .sub-header {
        font-size: 1.5rem;
        color: #5C7AEA;
        margin-bottom: 1rem;
    }
    .info-text {
        background-color: #f0f2f6;
        padding: 1rem;
        border-radius: 5px;
        margin-bottom: 1rem;
    }
    .success-box {
        background-color: #d4edda;
        color: #155724;
        padding: 1rem;
        border-radius: 5px;
        margin-bottom: 1rem;
    }
    .stProgress .st-eb {
        background-color: #4257B2;
    }
</style>
""", unsafe_allow_html=True)

# Header
st.markdown("<h1 class='main-header'>🧬 Protein Structure Prediction App</h1>", unsafe_allow_html=True)
st.markdown("<p>Predict and visualize protein structures using ESMFold</p>", unsafe_allow_html=True)

@st.cache_data
def render_mol(pdb_str):
    """Render molecule from PDB string"""
    viewer = py3Dmol.view(width=700, height=500)
    viewer.addModel(pdb_str, 'pdb')
    viewer.setStyle({'cartoon': {'color': 'spectrum'}})
    viewer.zoomTo()
    viewer.spin(True)
    viewer.zoom(0.8)
    return viewer

def render_from_pdb_file(pdb_file):
    """Render molecule from PDB file"""
    pdb_str = pdb_file.getvalue().decode('utf-8')
    return render_mol(pdb_str)

def predict_structure(sequence):
    """Predict protein structure using the ESM Metagenomic Atlas API"""
    url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
    headers = {"Content-Type": "application/x-www-form-urlencoded"}
    
    try:
        response = requests.post(url, headers=headers, data=sequence, timeout=30)
        if response.status_code == 200:
            return response.text
        else:
            st.error(f"Error: {response.status_code} - {response.text}")
            return None
    except requests.exceptions.RequestException as e:
        st.error(f"Request failed: {e}")
        return None

def calculate_protein_properties(sequence):
    """Calculate basic protein properties"""
    # Amino acid molecular weights
    aa_weights = {
        'A': 89.09, 'R': 174.20, 'N': 132.12, 'D': 133.10, 'C': 121.16,
        'Q': 146.15, 'E': 147.13, 'G': 75.07, 'H': 155.16, 'I': 131.17,
        'L': 131.17, 'K': 146.19, 'M': 149.21, 'F': 165.19, 'P': 115.13,
        'S': 105.09, 'T': 119.12, 'W': 204.23, 'Y': 181.19, 'V': 117.15
    }
    
    # Count amino acids
    aa_count = {}
    for aa in sequence:
        if aa in aa_count:
            aa_count[aa] += 1
        else:
            aa_count[aa] = 1
    
    # Calculate molecular weight
    mol_weight = sum(aa_weights.get(aa, 0) * count for aa, count in aa_count.items())
    
    # Calculate isoelectric point (simplified estimation)
    # This is a simplified calculation and not accurate for all proteins
    basic_aa = sequence.count('R') + sequence.count('K') + sequence.count('H')
    acidic_aa = sequence.count('D') + sequence.count('E')
    pi_est = 7.0 + 0.1 * (basic_aa - acidic_aa)
    pi_est = max(1.0, min(14.0, pi_est))  # Clamp to valid pH range
    
    return {
        "length": len(sequence),
        "molecular_weight": round(mol_weight, 2),
        "estimated_pi": round(pi_est, 2),
        "aa_composition": aa_count
    }

def show_protein_info(sequence, properties):
    """Display protein information"""
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("<h3 class='sub-header'>Protein Properties</h3>", unsafe_allow_html=True)
        st.info(f"Sequence Length: {properties['length']} amino acids")
        st.info(f"Molecular Weight: {properties['molecular_weight']} Da")
        st.info(f"Estimated Isoelectric Point: {properties['estimated_pi']}")
    
    with col2:
        st.markdown("<h3 class='sub-header'>Amino Acid Composition</h3>", unsafe_allow_html=True)
        
        # Count amino acids
        aa_count = properties['aa_composition']
        
        # Calculate percentages
        total_aa = sum(aa_count.values())
        aa_percentages = {aa: (count / total_aa) * 100 for aa, count in aa_count.items()}
        
        # Sort by count
        sorted_aa = dict(sorted(aa_percentages.items(), key=lambda x: x[1], reverse=True))
        
        # Display as bar chart
        st.bar_chart(sorted_aa)

def is_valid_sequence(sequence):
    """Check if the sequence contains only valid amino acid letters"""
    valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
    return all(aa in valid_aa for aa in sequence.upper())

def app():
    """Main application function"""
    tabs = st.tabs(["Predict Structure", "About ESMFold", "Help"])
    
    with tabs[0]:
        st.markdown("<h2 class='sub-header'>Input Protein Sequence</h2>", unsafe_allow_html=True)
        
        input_method = st.radio(
            "Choose input method:",
            ("Enter protein sequence", "Upload FASTA file", "Use example sequence")
        )
        
        sequence = ""
        
        if input_method == "Enter protein sequence":
            sequence = st.text_area(
                "Enter protein sequence (one-letter amino acid code):",
                height=150,
                help="Enter a protein sequence using the standard one-letter amino acid code (ACDEFGHIKLMNPQRSTVWY)"
            ).strip().upper().replace(" ", "").replace("\n", "")
            
        elif input_method == "Upload FASTA file":
            fasta_file = st.file_uploader("Upload a FASTA file", type=["fasta", "fa", "txt"])
            if fasta_file:
                fasta_content = fasta_file.getvalue().decode("utf-8")
                lines = fasta_content.strip().split("\n")
                # Skip header lines that start with >
                sequence_lines = [line.strip() for line in lines if not line.startswith(">")]
                sequence = "".join(sequence_lines).upper()
                
        elif input_method == "Use example sequence":
            example = st.selectbox(
                "Choose an example protein:",
                ("Short peptide (Bradykinin)", "Small protein (Insulin B chain)", "Medium protein (Lysozyme fragment)")
            )
            
            if example == "Short peptide (Bradykinin)":
                sequence = "RPPGFSPFR"
            elif example == "Small protein (Insulin B chain)":
                sequence = "FVNQHLCGSHLVEALYLVCGERGFFYTPKT"
            elif example == "Medium protein (Lysozyme fragment)":
                sequence = "KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL"
        
        # Display sequence info
        if sequence:
            st.markdown("<div class='info-text'>", unsafe_allow_html=True)
            st.write(f"Sequence length: {len(sequence)} amino acids")
            st.markdown("</div>", unsafe_allow_html=True)
            
            if not is_valid_sequence(sequence):
                st.error("Invalid sequence! Please use only standard amino acid letters (ACDEFGHIKLMNPQRSTVWY).")
                sequence = "" 
        
        # Button to predict structure
        if sequence:
            col1, col2 = st.columns([1, 3])
            
            with col1:
                predict_button = st.button("Predict Structure", type="primary", use_container_width=True)
            
            with col2:
                if len(sequence) > 400:
                    st.warning("⚠️ Warning: Long sequences (>400 amino acids) may take several minutes to complete or potentially fail.")
            
            if predict_button:
                if len(sequence) > 1000:
                    st.error("Error: Sequence is too long. Please use sequences with fewer than 1000 amino acids.")
                else:
                    with st.spinner("Predicting protein structure... This may take a few minutes."):
                        # Show progress bar
                        progress_bar = st.progress(0)
                        for i in range(100):
                            # Simulating work being done
                            time.sleep(0.05)
                            progress_bar.progress(i + 1)
                        
                        # Calculate properties
                        properties = calculate_protein_properties(sequence)
                        show_protein_info(sequence, properties)
                        
                        # Predict structure
                        pdb_str = predict_structure(sequence)
                        
                        if pdb_str:
                            # Show success message
                            st.markdown("<div class='success-box'>✅ Structure prediction successful!</div>", unsafe_allow_html=True)
                            
                            # Display the structure
                            st.markdown("<h3 class='sub-header'>Predicted 3D Structure</h3>", unsafe_allow_html=True)
                            view = render_mol(pdb_str)
                            showmol(view, height=500, width=800)
                            
                            # Download options
                            st.download_button(
                                label="Download PDB File",
                                data=pdb_str,
                                file_name="predicted_structure.pdb",
                                mime="chemical/x-pdb"
                            )
    
    with tabs[1]:
        st.markdown("<h2 class='sub-header'>About ESMFold</h2>", unsafe_allow_html=True)
        
        st.markdown("""
        ESMFold is a state-of-the-art protein structure prediction model developed by Meta AI Research. It's based on the ESM-2 language model and can accurately predict protein structures directly from amino acid sequences.
        
        **Key features of ESMFold:**
        
        * Predicts protein structures end-to-end using a language model approach
        * Fast inference compared to traditional molecular dynamics approaches
        * Competitive accuracy with AlphaFold2 for many protein structures
        * Can predict structures for sequences with no known homologs
        
        This app uses the ESM Metagenomic Atlas API to generate predictions.
        
        **References:**
        
        Lin, Z., Akin, H., Rao, R., et al. (2023). Evolutionary-scale prediction of atomic-level protein structure with a language model. Science, 379(6637), 1123-1130.
        """)
        
        col1, col2 = st.columns(2)
        with col1:
            st.image("https://raw.githubusercontent.com/facebookresearch/esm/main/assets/esmfold-plddt.gif", caption="ESMFold prediction confidence map (source: Meta AI)", use_column_width=True)
        
    with tabs[2]:
        st.markdown("<h2 class='sub-header'>Help & Instructions</h2>", unsafe_allow_html=True)
        
        st.markdown("""
        ### How to use this app:
        
        1. **Input your protein sequence** using one of three methods:
           - Type or paste a protein sequence in one-letter amino acid code
           - Upload a FASTA file containing the sequence
           - Select one of the example sequences
           
        2. **Click "Predict Structure"** to start the prediction process
        
        3. **View the results**:
           - The 3D structure visualization will appear after successful prediction
           - You can rotate, zoom, and explore the structure interactively
           - Download the PDB file for use in other molecular visualization software
           
        ### Tips for best results:
        
        - For optimal performance, use sequences between 50-400 amino acids
        - Very short sequences (<20 aa) may not fold into stable structures
        - Very long sequences (>400 aa) may take longer to process
        - Make sure your sequence contains only standard amino acid letters
        
        ### Troubleshooting:
        
        - If prediction fails, try again or use a shorter sequence
        - Check that your sequence contains only valid amino acid letters
        - For large proteins, consider splitting into domains or using other tools
        """)

if __name__ == "__main__":
    app()
