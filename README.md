# M.A.P.P.E.R. 🧬
**Motif Analyzer and Promoter Predictor for Expressions and Regulation**

M.A.P.P.E.R. is a local web application developed in Python (using Streamlit) designed for the advanced analysis of bacterial genomes. It automates the discovery of Open Reading Frames (ORFs) and predicts the complete structure of operons by combining the detection of promoter sequences (Sigma factors), transcriptional regulator binding sites, Ribosome Binding Sites (RBS / Shine-Dalgarno), and intrinsic terminators.

Initially calibrated with default Position Weight Matrices (PWM/PSSM) for the *Bacillus* genus, the tool is universal and supports custom matrix imports to scan any bacterial contig from sequencing data (FASTA/Multi-FASTA formats).

## 🚀 Features
* **Multi-FASTA Support**: Full automated analysis across all contigs or fragments within a single sequencing file.
* **Linear-Time Engine $O(N)$**: Optimized algorithm capable of scanning large genomic fragments (mega-plasmids, chromosomes) without memory overhead.
* **Real-Time Dynamic Filtering**: Instant adjustment of affinity thresholds via UI sliders without triggering a full genomic recalculation (powered by a dual-stage extraction architecture).
* **Multi-Fragment Navigation**: Intuitive dropdown menu to seamlessly switch between different contigs within the same Multi-FASTA file.
* **Confidence Scoring**: Dynamic scoring system (/100) for each gene and operon based on biological signals (quality of -35/-10 boxes, RBS affinity, GC3 bias, and terminator thermostability).
* **Zoomed Visual Mapping**: Clean, BioRender-style vector genomic maps centered on detected operons with flanging margins to assist in cloning primer design.
* **Alternative Start Verification**: Dropdown menus to inspect alternative initiation codons upstream/downstream, ensuring the safety of native promoter amplifications.
* **Custom Matrices**: Upload organism-specific regulator sequences, Sigma factors, and custom spacer ranges using a standard Excel `.xlsx` file.
* **Excel Exporting**: One-click export of complete data tables (including coordinates, scores, annotations, and translated protein sequences) for future experimental planning.

## 🛠️ Installation & Prerequisites

The application runs entirely locally on your machine. No internet connection is required after installation, guaranteeing the strict confidentiality of your genomic data. You will need **Python 3.8+** installed.

1. **Clone the repository:**
  `git clone [https://github.com/your-name/MAPPER.git](https://github.com/your-name/MAPPER.git)`
  `cd MAPPER`

**Install the required dependencies:**
It is highly recommended to use a virtual environment (e.g., via Conda or venv). Install the dependencies using:
`pip install -r requirements.txt`

## 🖥️ Usage
To launch the application, open your terminal in the project directory and run:
`streamlit run sequence_analyzer.py`
A window will automatically open in your default web browser at the local address `http://localhost:8501`

**Analysis Steps:**
1. Configuration: Adjust the minimum ORF length and initial affinity thresholds for promoters, RBS, and regulators in the sidebar.

2. Matrices (Optional): If working on a different biological model, download the pre-filled Excel template from the UI, edit it with your motifs of interest, and re-upload it. If left blank, the default Bacillus dictionnaire will be automatically applied.

3. Scan: Drag and drop your genomic file (.fasta, .fa, .fna). The initial run will extract all potential features.

4. Explore & Fine-Tune: Use the sliders to filter out background noise in real-time. Navigate through your contigs and use the pagination to visually inspect operon maps.

5. Export: Click "Download full table (Excel)" to retrieve the final summary table containing all annotations and sequences for your downstream benchwork.

6. Export : Une fois le calcul terminé, naviguez visuellement contig par contig ou cliquez sur "Télécharger le tableau complet (Excel)" pour extraire toutes vos données.
