import streamlit as st  # type: ignore
from Bio import SeqIO  # type: ignore
import pandas as pd  # type: ignore
import io
import plotly.express as px  # type: ignore
import py3Dmol  # type: ignore
import tempfile
import os

# Title of the app
st.title("¡Sube tu archivo FASTA aquí!")
st.write("Por favor, sube tu archivo FASTA para comenzar el análisis.")

# File uploader: The user will upload their FASTA file here
uploaded_fasta_file = st.file_uploader("Elige tu archivo FASTA", type=["fasta", "fa"])

if uploaded_fasta_file is not None:
    # Convert uploaded file to a text stream
    stringio = io.StringIO(uploaded_fasta_file.getvalue().decode("utf-8"))
    
    # Parse the FASTA file from the text stream
    record = SeqIO.read(stringio, "fasta")
    
    # Show sequence ID and sequence as text
    st.subheader(f"ID de la secuencia: {record.id}")
    st.text_area("Secuencia", str(record.seq), height=200)

    # Selection for visualization methods
    st.subheader("Selecciona método(s) de visualización")
    show_3d_model = st.checkbox("Modelo 3D")
    show_table = st.checkbox("Tabla de secuencia")
    show_aa_distribution = st.checkbox("Distribución de Aminoácidos")
    show_sequence_alignment = st.checkbox("Alineamiento de Secuencias")
    show_secondary_structure = st.checkbox("Estructura Secundaria de la Proteína")
    
    # Display methods based on checkbox selection
    if show_table:
        # Count the number of repetitions for each amino acid in the sequence
        aa_count = {aa: str(record.seq).count(aa) for aa in set(str(record.seq))}
        aa_df = pd.DataFrame(list(aa_count.items()), columns=["Aminoácido", "Repeticiones"])
        aa_df = aa_df.sort_values(by="Repeticiones", ascending=False)  # Sort by count
        st.write("Tabla de Aminoácidos y sus Repeticiones")
        st.write(aa_df)

    if show_aa_distribution:
        # Count amino acids and their percentage
        aa_count = {aa: str(record.seq).count(aa) for aa in set(str(record.seq))}
        aa_labels = list(aa_count.keys())
        aa_values = list(aa_count.values())
        total_aa = sum(aa_values)
        aa_percentages = [round((value / total_aa) * 100, 2) for value in aa_values]

        # Create a DataFrame for Plotly
        aa_df = pd.DataFrame({
            "Aminoácido": aa_labels,
            "Repeticiones": aa_values,
            "Porcentaje": aa_percentages
        })

        # Create a pie chart with Plotly
        fig = px.pie(aa_df, names='Aminoácido', values='Porcentaje', hover_data={'Aminoácido': True, 'Porcentaje': True},
                     title="Distribución de Aminoácidos", 
                     labels={"Porcentaje": "Porcentaje (%)", "Aminoácido": "Aminoácido"})
        
        # Customize hover info to display the name and percentage
        fig.update_traces(textinfo='percent+label', hovertemplate='%{label}: %{value}%')
        
        # Show the pie chart
        st.subheader("Distribución de Aminoácidos")
        st.plotly_chart(fig)

    if show_sequence_alignment:
        st.write("Puedes realizar un alineamiento de secuencias usando herramientas como [Clustal Omega](https://www.ebi.ac.uk/Tools/msa/clustalo/) o [MAFFT](https://mafft.cbrc.jp/alignment/software/).")
        
        # Create a file uploader specifically for sequence alignment
        alignment_file = st.file_uploader("Sube tu archivo para el alineamiento de secuencias (FASTA)", type=["fasta", "fa"])
        if alignment_file is not None:
            alignment_stringio = io.StringIO(alignment_file.getvalue().decode("utf-8"))
            alignment_record = SeqIO.read(alignment_stringio, "fasta")
            st.write(f"Secuencia cargada para el alineamiento: {alignment_record.id}")
            st.text_area("Secuencia cargada", str(alignment_record.seq), height=200)

    if show_secondary_structure:
        st.write("Puedes predecir la estructura secundaria de la proteína usando herramientas como [PSIPRED](https://www.ebi.ac.uk/Tools/psipred/).")
        
        # Create a file uploader for secondary structure prediction
        secondary_structure_file = st.file_uploader("Sube tu archivo para la predicción de la estructura secundaria (FASTA)", type=["fasta", "fa"])
        if secondary_structure_file is not None:
            # Procesar el archivo cargado para predicción
            secondary_structure_stringio = io.StringIO(secondary_structure_file.getvalue().decode("utf-8"))
            secondary_structure_record = SeqIO.read(secondary_structure_stringio, "fasta")
            st.write(f"Secuencia cargada para predicción de estructura secundaria: {secondary_structure_record.id}")
            st.text_area("Secuencia cargada", str(secondary_structure_record.seq), height=200)

    if show_3d_model:
        # Sección para cargar un archivo PDB y visualizar el modelo 3D
        st.write("Para generar un modelo 3D a partir de tu secuencia de aminoácidos, puedes utilizar las siguientes herramientas:")
        st.markdown("[*I-TASSER*](https://zhanggroup.org/I-TASSER/) - Para generar un modelo 3D a partir de tu secuencia de proteínas.")
        st.markdown("[*AlphaFold*](https://alphafold.ebi.ac.uk/) - Predicción de la estructura 3D usando inteligencia artificial de DeepMind.")
        st.write("Una vez generado el modelo 3D, puedes cargar el archivo PDB aquí para visualizarlo.")
        
        # Create a file uploader specifically for the PDB file to display 3D model
        pdb_file = st.file_uploader("Sube tu archivo PDB para visualizar el modelo 3D", type=["pdb"])
        
        if pdb_file is not None:
            try:
                # Save the uploaded PDB file to a temporary location
                with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as tmp_file:
                    tmp_file.write(pdb_file.getvalue())
                    tmp_file_path = tmp_file.name  # Temporary file path
                
                # Ensure the file is a valid PDB file
                if tmp_file_path.endswith('.pdb'):
                    # Read the PDB file
                    with open(tmp_file_path, 'r') as file:
                        pdb_data = file.read()
                    
                    # Visualize using py3Dmol
                    viewer = py3Dmol.view(width=800, height=400)
                    viewer.addModel(pdb_data, "pdb")
                    viewer.setStyle({'cartoon': {'color': 'spectrum'}})  # Style the structure
                    viewer.zoomTo()
                    
                    # Display the 3D model in Streamlit
                    viewer.show()
                    
                else:
                    st.error("Por favor, sube un archivo con extensión .pdb")
                
                # Clean up temporary file
                os.remove(tmp_file_path)
            except Exception as e:
                st.write(f"Error al cargar el archivo PDB: {e}")

