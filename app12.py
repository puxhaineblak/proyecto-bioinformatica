import streamlit as st
from Bio import SeqIO
import pandas as pd
import io
import plotly.express as px
import py3Dmol
import tempfile
import os

# Title of the app
st.title("¡Sube tu archivo FASTA aquí!")
st.write("Por favor, sube tu archivo FASTA para comenzar el análisis.")

# File uploader: The user will upload their FASTA file here
uploaded_file = st.file_uploader("Elige tu archivo FASTA", type=["fasta", "fa"])

def count_amino_acids(sequence):
    aa_count = {aa: sequence.count(aa) for aa in set(sequence)}
    return aa_count

if uploaded_file is not None:
    # Convert uploaded file to a text stream and parse the FASTA file
    try:
        with io.StringIO(uploaded_file.getvalue().decode("utf-8")) as stringio:
            record = SeqIO.read(stringio, "fasta")
    except Exception as e:
        st.error(f"Error al leer el archivo FASTA: {e}")
        record = None

    if record:
        st.subheader(f"ID de la secuencia: {record.id}")
        st.text_area("Secuencia", str(record.seq), height=200)

        # Selection for visualization methods
        st.subheader("Selecciona método(s) de visualización")
        show_3d_model = st.checkbox("Modelo 3D")
        show_table = st.checkbox("Tabla de secuencia")
        show_aa_distribution = st.checkbox("Distribución de Aminoácidos")
        show_sequence_alignment = st.checkbox("Alineamiento de Secuencias")
        show_secondary_structure = st.checkbox("Estructura Secundaria de la Proteína")

        if show_table:
            aa_count = count_amino_acids(str(record.seq))
            aa_df = pd.DataFrame(list(aa_count.items()), columns=["Aminoácido", "Repeticiones"])
            aa_df = aa_df.sort_values(by="Repeticiones", ascending=False)
            st.write("Tabla de Aminoácidos y sus Repeticiones")
            st.write(aa_df)

        if show_aa_distribution:
            aa_count = count_amino_acids(str(record.seq))
            aa_labels = list(aa_count.keys())
            aa_values = list(aa_count.values())
            total_aa = sum(aa_values)
            aa_percentages = [round((value / total_aa) * 100, 2) for value in aa_values]

            aa_df = pd.DataFrame({
                "Aminoácido": aa_labels,
                "Repeticiones": aa_values,
                "Porcentaje": aa_percentages
            })

            fig = px.pie(aa_df, names='Aminoácido', values='Porcentaje', hover_data={'Aminoácido': True, 'Porcentaje': True},
                         title="Distribución de Aminoácidos", 
                         labels={"Porcentaje": "Porcentaje (%)", "Aminoácido": "Aminoácido"})
            fig.update_traces(textinfo='percent+label', hovertemplate='%{label}: %{value}%')
            st.subheader("Distribución de Aminoácidos")
            st.plotly_chart(fig)

        if show_sequence_alignment:
            st.write("Puedes realizar un alineamiento de secuencias usando herramientas como [Clustal Omega](https://www.ebi.ac.uk/Tools/msa/clustalo/) o [MAFFT](https://mafft.cbrc.jp/alignment/software/)")
            alignment_file = st.file_uploader("Sube tu archivo para el alineamiento de secuencias (FASTA)", type=["fasta", "fa"])
            if alignment_file is not None:
                try:
                    with io.StringIO(alignment_file.getvalue().decode("utf-8")) as alignment_stringio:
                        alignment_record = SeqIO.read(alignment_stringio, "fasta")
                    st.write(f"Secuencia cargada para el alineamiento: {alignment_record.id}")
                    st.text_area("Secuencia cargada", str(alignment_record.seq), height=200)
                except Exception as e:
                    st.error(f"Error al leer el archivo de alineamiento: {e}")

        if show_secondary_structure:
            st.write("Puedes predecir la estructura secundaria de la proteína usando herramientas como [PSIPRED](https://www.ebi.ac.uk/Tools/psipred/).")
            secondary_structure_file = st.file_uploader("Sube tu archivo para la predicción de la estructura secundaria (FASTA)", type=["fasta", "fa"])
            if secondary_structure_file is not None:
                try:
                    with io.StringIO(secondary_structure_file.getvalue().decode("utf-8")) as secondary_structure_stringio:
                        secondary_structure_record = SeqIO.read(secondary_structure_stringio, "fasta")
                    st.write(f"Secuencia cargada para predicción de estructura secundaria: {secondary_structure_record.id}")
                    st.text_area("Secuencia cargada", str(secondary_structure_record.seq), height=200)
                except Exception as e:
                    st.error(f"Error al leer el archivo de estructura secundaria: {e}")

        if show_3d_model:
            st.write("Para generar un modelo 3D a partir de tu secuencia de aminoácidos, puedes utilizar las siguientes herramientas:")
            st.markdown("[*I-TASSER*](https://zhanggroup.org/I-TASSER/) - Para generar un modelo 3D a partir de tu secuencia de proteínas.")
            st.markdown("[*AlphaFold*](https://alphafold.ebi.ac.uk/) - Predicción de la estructura 3D usando inteligencia artificial de DeepMind.")
            st.write("Una vez generado el modelo 3D, puedes cargar el archivo PDB aquí para visualizarlo.")
            
            pdb_file = st.file_uploader("Sube tu archivo PDB para visualizar el modelo 3D", type=["pdb"])
            if pdb_file is not None:
                try:
                    with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as tmp_file:
                        tmp_file.write(pdb_file.getvalue())
                        tmp_file_path = tmp_file.name
                    
                    if tmp_file_path.endswith('.pdb'):
                        with open(tmp_file_path, 'r') as file:
                            pdb_data = file.read()
                        viewer = py3Dmol.view(width=800, height=400)
                        viewer.addModel(pdb_data, "pdb")
                        viewer.setStyle({'cartoon': {'color': 'spectrum'}})
                        viewer.zoomTo()
                        viewer.show()
                    else:
                        st.error("Por favor, sube un archivo con extensión .pdb")
                except Exception as e:
                    st.write(f"Error al cargar el archivo PDB: {e}")
                finally:
                    if os.path.exists(tmp_file_path):
                        os.remove(tmp_file_path)

