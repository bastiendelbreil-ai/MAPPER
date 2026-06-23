import re
import copy
from io import BytesIO, StringIO
from itertools import zip_longest

import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
from Bio import SeqIO, motifs
from Bio.Seq import Seq
from dna_features_viewer import GraphicFeature, GraphicRecord


# ==========================================
# 1. GLOBAL CONFIGURATION & PARAMETERS
# ==========================================

class Config:
    """Centralizes all distances and biological parameters."""
    MAX_PROMOTER_DIST = 500
    RBS_ZONE = 30
    TERM_ZONE = 100
    PROM_ZONE = 300
    MAX_GAP = 50
    MAX_OVERLAP = 40
    START_CODONS = ["ATG", "TTG", "GTG"]
    STOP_CODONS = ["TAA", "TAG", "TGA"]

# Pre-compiled Regex for faster execution
DIGIT_PATTERN = re.compile(r'\d+')

DEFAULT_REGULATORS = {
    "CodY": [
        Seq("AATTTTCAGAAAATT"), Seq("AATTTTCTGAAAATT"), Seq("AACTTTCAGAATATT"), Seq("AATATTCTAATAATT"), Seq("GATTTTCTGAAAATC"), Seq("ATTTTTCAGAAAATC"), 
        Seq("ATTTTTCTAAATAGT"), Seq("ATTTTCCTAAATAGT"), Seq("ATTTTTCTGAAAATT"), Seq("ATTTTTCTAATAATT")
    ],
    "PlcR": [
        Seq("TATGAAACAATTCATA"), Seq("TATGGAAAGGTCCATA"), Seq("TATGCAATTATACATA"),
        Seq("TATGCAAAAATGCATA"), Seq("TATGAGAAAATTCATA")  
    ],
    "CRE high affinity": [
        Seq("TGAAAGCGCTTTAT"), Seq("TGAAAGCGTTACCA"), Seq("TGAAAACGCTTTAT"), Seq("TGTAAGCGTTAACA"), 
        Seq("TGAAAGCGTTGACA"), Seq("TGTAAGCGTATACA"), Seq("TGTAAGCGGATTCA"), Seq("TGAATGCGGTTACA"), Seq("TGAAAACGCTATCA"), 
        Seq("TGACACCGCTTTCA"), Seq("TGTAAGCGTTTTAA"), Seq("TGAAAGCGTTTAAT"), Seq("TGTAAACGGTTTCT"), Seq("AGAAAGCGTTTACA"),  
        Seq("TGTAAACGGTTACA"), Seq("TTAAAGCGCTTACA"), Seq("CGAAAACGCTATCA"), Seq("TGAAAGCGCAGTCT"), Seq("TGAAAACGCTTGCA"),         
        Seq("TGAAAGCGTTATCA"), Seq("TAAAAGCGCTTACA"), Seq("TGGAAGCGCAAACA"), Seq("TGAAAGCGCTTTTA"), Seq("AGAAAACGCTTTCT"),       
        Seq("TGAAAGCGGTTCAA"), Seq("TGCAAGGGCTTTCA"), Seq("TGATAGCGTTTTCT"), Seq("TGTAAGCGTTTGCT"), Seq("TGAAAGCGCTATCT")        
    ],
    "CRE low affinity": [
        Seq("TGTAAGCGCTTTCT"), Seq("TGTAAGCGTTTGCT"), Seq("TGTAAGCATTTTCT"), Seq("TGAAAACGGTGCCA"), Seq("TGAAATGAATCGTT"), 
        Seq("TGAAAACGGTTTAT"), Seq("TGAAAGTGTTTGCA"), Seq("TGAAAGCGGTACCA"), Seq("TGAAACCGCTTCCA"), Seq("TGAAAACGTTATCA"), 
        Seq("TGAAAACGCTTAAC"), Seq("TGTAAACGTTATCA"), Seq("TGAAAGCGTTTTAG"), Seq("TGTAACCGCTTACT"), Seq("TGAAAGCGTTTTAT"),  
        Seq("TGGAAGCGTTTTTA"), Seq("TGAAAACGTTTTTT"), Seq("TGAAAGCGCTATAA"), Seq("TAAAAACGCTTTCT"), Seq("GGAAAACGCTTTCA"),         
        Seq("TGTAAGCGGTATCT"), Seq("TGAAAACGCGATCA"), Seq("AGAAAGCGCTTACG"), Seq("TGAAAACGTGTCCA"), Seq("TGAAAGCGTTTTCC"),       
        Seq("AGAAAGCGGTTACA"), Seq("TGTAAACGGTTACA"), Seq("TGAAAGCGCTTTCT"), Seq("TGAAAACGCTATCG"), Seq("TGAAAATGTTTACA")          
    ],
    "CRE activator": [
        Seq("TGAAAGCGTATACA"), Seq("TGAAAGCGTTTTAT"), Seq("TGAAAGCGCTCGCT")          
    ],
    "Spo0A": [
        Seq("TGTCGAA")   
    ],
    "RBS (Shine-Dalgarno)": [
        Seq("AGGAGG"), Seq("GGAGGG"), Seq("AAGGAG"), Seq("GAGGGG"), Seq("GAAAGG"), Seq("AAGGGG"), Seq("GGGTGG"), Seq("GGGGGA")
    ],
    "SigmaA_-35": [
        Seq("TGTACA"), Seq("TAGTAA"), Seq("TTGACG"), Seq("TTGAAA"), Seq("TTGAGA"), Seq("TTGTTA"), Seq("TTGCAA"),
        Seq("TGGTCA"), Seq("TTGGAG"), Seq("TTGCAG"), Seq("TGGACA"), Seq("TGGAGG"), Seq("CTGTCG"), Seq("TTGTTT"),
        Seq("TTGCGT"), Seq("TTGTAT"), Seq("TTGTCC"), Seq("TTGCCA"), Seq("TCGAAA"), Seq("CTGACT"), Seq("TTGTCA"),
        Seq("TTGAGT"), Seq("TTGCAT"), Seq("TGGTAT"), Seq("TAGACA"), Seq("TTGACA"), Seq("TTGACC"), Seq("TTGGCC"),
        Seq("TTGACT"), Seq("TTGTCT"), Seq("TGGTGT"), Seq("TTGGAA"), Seq("GTGAAA"), Seq("TTGAAT"), Seq("TTGATT")
    ],
    "SigmaA_-10": [
        Seq("TACAAT"), Seq("TATGCT"), Seq("TATAAT"), Seq("TATTCT"), Seq("TATAGT"), Seq("TAAAAT"), Seq("TACACT"),
        Seq("TAAACT"), Seq("TATATT"), Seq("TAATTT"), Seq("TAGAGT"), Seq("TATACT"), Seq("TATGAT"), Seq("TACAGT"),
        Seq("TACTAT"), Seq("TATCAT"), Seq("TACCAT"), Seq("TATAAA"), Seq("CATACT"), Seq("TAGATT"), Seq("TATTAT"),
        Seq("TAGAAT"), Seq("TAAAGT"), Seq("TAGACT"), Seq("TATACA"), Seq("TACTCT")
    ],
   "SigmaB_-35": [
        Seq("GTTTAA"), Seq("GTTTTA"), Seq("GTTTCA"), Seq("GTTTGA"), Seq("GTTTAT")
    ],
    "SigmaB_-10": [
        Seq("GGGTAT"), Seq("GGGTAC"), Seq("GGGTTA"), Seq("GGGTGT"), Seq("GGATAT")
    ],
    "SigmaE_-35": [
        Seq("GCATAT"), Seq("GCATAA"), Seq("GCATAC"), Seq("ATATAT"), 
        Seq("GCAAAT"), Seq("GCATGT"), Seq("GCATTA"), Seq("TCATAT")
    ],
    "SigmaE_-10": [
        Seq("CATACA"), Seq("CATACT"), Seq("CATAAT"), Seq("TATACA"), 
        Seq("CATACC"), Seq("CATAAA")
    ],
    "SigmaK_-35": [
        Seq("GCAGTA"), Seq("GCAGCA"), Seq("CCAGTA"), Seq("TCAGTA"), 
        Seq("GCAATA"), Seq("GCAGAA"), Seq("GCAGTC"), Seq("GCACGA")
    ],
    "SigmaK_-10": [
        Seq("CATGTA"), Seq("CATAAT"), Seq("CATACA"), Seq("CATAAA"), 
        Seq("TATGTA"), Seq("CATGAA")
    ],
    "SigmaG_-35": [
        Seq("GCAGGA"), Seq("GCAAGA"), Seq("GCACGA"), Seq("CCAGGA"), 
        Seq("TCAGGA"), Seq("GCAGTA")
    ],
    "SigmaG_-10": [
        Seq("CATACA"), Seq("CATAAT"), Seq("CATACT"), Seq("TATACA"), 
        Seq("CATAAA"), Seq("CATGCA")
    ]
}

DEFAULT_SPACERS = {
    "SigmaA": [15, 16, 17, 18],
    "SigmaB": [13, 14, 15],
    "SigmaE": [13, 14, 15],
    "SigmaK": [13, 14, 15],
    "SigmaG": [15, 16, 17, 18]
}

# ==========================================
# 2. CORE BIOINFORMATICS FUNCTIONS
# ==========================================

def decode_spacer(cell_value):
    if pd.isna(cell_value) or str(cell_value).strip() == "":
        return [15, 16, 17, 18]
        
    val_str = str(cell_value).strip()
    numbers = [int(n) for n in DIGIT_PATTERN.findall(val_str)]
    
    if len(numbers) == 1:
        return [numbers[0]]
    elif len(numbers) >= 2:
        return list(range(min(numbers), max(numbers) + 1))
    return [15, 16, 17, 18]

def read_custom_matrices(excel_file):
    custom_matrices = {}
    custom_spacers = {}
    
    try:
        df_prom = pd.read_excel(excel_file, sheet_name="Promoters")
        for _, row in df_prom.iterrows():
            sigma_name = str(row.iloc[0]).strip()
            seq_35 = str(row.iloc[1]).strip().upper()
            seq_10 = str(row.iloc[2]).strip().upper()
            
            if pd.isna(sigma_name) or sigma_name == "nan":
                continue
            
            try:
                spacer_val = row.iloc[3]
            except IndexError:
                spacer_val = None
                
            custom_spacers[sigma_name] = decode_spacer(spacer_val)
            key_35 = f"{sigma_name}_-35"
            key_10 = f"{sigma_name}_-10"
            
            if key_35 not in custom_matrices: custom_matrices[key_35] = []
            if key_10 not in custom_matrices: custom_matrices[key_10] = []
            
            if len(seq_35) >= 3 and seq_35 != "NAN":
                custom_matrices[key_35].append(Seq(seq_35))
            if len(seq_10) >= 3 and seq_10 != "NAN":
                custom_matrices[key_10].append(Seq(seq_10))
    except Exception:
        st.error("❌ Failed to read 'Promoters' sheet.")

    try:
        df_reg = pd.read_excel(excel_file, sheet_name="Regulators")
        for _, row in df_reg.iterrows():
            reg_name = str(row.iloc[0]).strip()
            seq_reg = str(row.iloc[1]).strip().upper()
            
            if pd.isna(reg_name) or reg_name == "nan" or len(seq_reg) < 3 or seq_reg == "NAN":
                continue
            if reg_name not in custom_matrices:
                custom_matrices[reg_name] = []
            custom_matrices[reg_name].append(Seq(seq_reg))
    except Exception:
        st.error("❌ Failed to read 'Regulators' sheet.")

    try:
        df_rbs = pd.read_excel(excel_file, sheet_name="RBS")
        custom_matrices["RBS (Shine-Dalgarno)"] = []
        for _, row in df_rbs.iterrows():
            seq_rbs = str(row.iloc[0]).strip().upper()
            
            if pd.isna(seq_rbs) or seq_rbs == "nan" or len(seq_rbs) < 3 or seq_rbs == "NAN":
                continue
            custom_matrices["RBS (Shine-Dalgarno)"].append(Seq(seq_rbs))
    except Exception:
        st.error("❌ Failed to read 'RBS' sheet.")
        
    return custom_matrices, custom_spacers

@st.cache_data
def build_pssm_matrices(_sequence_dict):
    pssms = {}
    for reg_name, sequences in _sequence_dict.items():
        try:
            m = motifs.create(sequences)
            pssm = m.counts.normalize(pseudocounts=0.5).log_odds()
            pssms[reg_name] = pssm
        except ValueError:
            st.warning(f"Warning: Cannot build matrix for '{reg_name}' (unequal lengths).")
    return pssms
    
def calculate_gc3(orf_seq):
    if len(orf_seq) < 3:
        return 0.0
    third_bases = orf_seq[2::3] 
    gc_count = third_bases.count('G') + third_bases.count('C')
    return (gc_count / len(third_bases)) * 100

def scan_regulators(promoter_seq, precalc_pssms, threshold_pct=0.70):
    found_sites = []
    seq_prom = Seq(promoter_seq)
    prom_len = len(seq_prom)
    
    for reg_name, pssm in precalc_pssms.items():
        max_score = pssm.max
        adapted_thresh = max_score * threshold_pct
        
        for pos, score in pssm.search(seq_prom, threshold=adapted_thresh):
            if pos >= 0: 
                rel_pos = pos - prom_len
                found_seq = seq_prom[pos : pos + pssm.length]
                real_pct = (score / max_score) * 100
                
                found_sites.append({
                    "Regulator": reg_name, "Sequence": str(found_seq),
                    "Position": f"{rel_pos} bp", "Score": round(score, 2),
                    "Percentage": round(real_pct, 1)
                })
    return found_sites

def find_rbs(upstream_seq, precalc_pssms, rbs_thresh):
    critical_zone = upstream_seq[-Config.RBS_ZONE:]
    zone_len = len(critical_zone)
    pssm_rbs = precalc_pssms["RBS (Shine-Dalgarno)"]
    max_score = pssm_rbs.max
    adapted_thresh = max_score * rbs_thresh
    
    found_rbs = [] 
    for pos, score in pssm_rbs.search(Seq(critical_zone), threshold=adapted_thresh):
        if pos >= 0:
            dist_to_start = zone_len - pos
            found_seq = critical_zone[pos : pos + pssm_rbs.length]
            pct = (score / max_score) * 100
            found_rbs.append({
                "Sequence": str(found_seq), "Distance": f"-{dist_to_start} bp",
                "Percentage": round(pct, 1), "Score": score
            })
            
    return sorted(found_rbs, key=lambda x: x['Score'], reverse=True)

def find_promoters(upstream_seq, precalc_pssms, prom_thresh, spacer_dict):
    prom_zone = upstream_seq[-Config.PROM_ZONE:]
    zone_len = len(prom_zone)
    found_promoters = []
    
    sigma_types = set(k.split('_')[0] for k in precalc_pssms.keys() if k.startswith("Sigma"))
    
    for sig_name in sigma_types:
        key_35 = f"{sig_name}_-35"
        key_10 = f"{sig_name}_-10"
        
        if key_35 not in precalc_pssms or key_10 not in precalc_pssms:
            continue
            
        pssm_35 = precalc_pssms[key_35]
        pssm_10 = precalc_pssms[key_10]
        max_total_score = pssm_35.max + pssm_10.max
        spacers = spacer_dict.get(sig_name, [15, 16, 17, 18])
        
        for pos_35, score_35 in pssm_35.search(Seq(prom_zone), threshold=0.0):
            if pos_35 >= 0:
                spacer_start = pos_35 + pssm_35.length
                for sp_len in spacers:
                    expected_10_pos = spacer_start + sp_len
                    if expected_10_pos + pssm_10.length <= zone_len:
                        found_10_seq = prom_zone[expected_10_pos : expected_10_pos + pssm_10.length]
                        score_10 = pssm_10.calculate(Seq(found_10_seq))
                        
                        total_score = score_35 + score_10
                        pct = (total_score / max_total_score) * 100
                        
                        if pct >= (prom_thresh * 100):
                            labels = {"SigmaA": "σA", "SigmaB": "σB", "SigmaE": "σE", "SigmaK": "σK", "SigmaG": "σG"}
                            found_promoters.append({
                                "Type": labels.get(sig_name, sig_name),
                                "Seq_35": str(prom_zone[pos_35 : pos_35 + 6]),
                                "Spacer": sp_len, "Seq_10": str(found_10_seq),
                                "Distance": f"-{zone_len - pos_35} bp",
                                "Percentage": round(pct, 1)
                            })
                            
    return sorted(found_promoters, key=lambda x: x['Percentage'], reverse=True)

def find_terminator(downstream_seq):
    best_term = None
    max_score = 0
    
    for i in range(len(downstream_seq) - 15):
        for stem_len in range(5, 13):
            for loop_len in range(3, 9):
                end_stem1 = i + stem_len
                end_loop = end_stem1 + loop_len
                end_stem2 = end_loop + stem_len
                
                if end_stem2 + 6 > len(downstream_seq):
                    continue
                    
                stem1 = downstream_seq[i:end_stem1]
                stem2 = downstream_seq[end_loop:end_stem2]
                tail = downstream_seq[end_stem2:end_stem2+6]
                
                stem1_rc = str(Seq(stem1).reverse_complement())
                matches = sum(1 for a, b in zip(stem1_rc, stem2) if a == b)
                
                if matches >= stem_len - 1:
                    t_count = tail.count('T')
                    if t_count >= 3:
                        gc_count = stem1.count('G') + stem1.count('C')
                        score = (matches * 2) + gc_count + (t_count * 1.5)
                        
                        if score > max_score:
                            max_score = score
                            best_term = {
                                "Stem1": stem1, "Loop": downstream_seq[end_stem1:end_loop],
                                "Stem2": stem2, "Tail": tail, "Distance": f"+{i} bp",
                                "Thermo_Score": round(score, 1)
                            }
    return best_term
    
def assign_confidence_score(orf, max_len):
    score = 0
    if len(orf['RBS']) > 0:
        best_rbs = orf['RBS'][0]
        rbs_pts = (best_rbs['Percentage'] / 100) * 25
        real_dist = int(best_rbs['Distance'].replace("-", "").replace(" bp", "").strip())
        if real_dist < 5 or real_dist > 15:
            rbs_pts *= 0.5 
        score += rbs_pts
            
    score += (orf['Length'] / max_len) * 25
            
    if orf['Terminator']:
        score += min(orf['Terminator']['Thermo_Score'], 15)
            
    if len(orf['Regulators']) > 0:
        sum_pts = sum(r['Percentage'] + (10 if r['Regulator'] == "CRE high affinity" else 0) for r in orf['Regulators'])
        avg_affinity = min(sum_pts / len(orf['Regulators']), 100.0) 
        score += (avg_affinity / 100) * 10
        
    if len(orf['Promoters']) > 0:
        score += (orf['Promoters'][0]['Percentage'] / 100) * 25
        
    if orf.get('GC3', 35) < 20 or orf.get('GC3', 35) > 55:
        score *= 0.80 
        
    score = min(score, 100)
    
    if score >= 80:   return "★★★★★", score
    elif score >= 65: return "★★★★☆", score
    elif score >= 50: return "★★★☆☆", score
    elif score >= 30: return "★★☆☆☆", score
    return "★☆☆☆☆", score

# ==========================================
# 3. MAIN ENGINE
# ==========================================

@st.cache_data
def extract_raw_orfs(seq_text, min_orf_len, _custom_dict, _custom_spacers):
    """ÉTAGE 1 : L'extraction lourde (Ne tourne qu'une seule fois)"""
    try:
        records = list(SeqIO.parse(StringIO(seq_text), "fasta"))
    except Exception as e:
        st.error(f"Error reading sequence file: {e}")
        return None, 0
        
    pssms = build_pssm_matrices(_custom_dict)
    all_raw_orfs = []
    total_analyzed_len = 0
    
    for idx_record, record in enumerate(records, 1):
        seq_obj = record.seq 
        
        nom_brut = record.description if record.description else record.id
        contig_id = f"{nom_brut} [Seq {idx_record}]"
        
        seq_len = len(seq_obj)
        total_analyzed_len += seq_len
        
        strands = [("Direct (+)", str(seq_obj).upper()), ("Reverse (-)", str(seq_obj.reverse_complement()).upper())]
        
        for strand_name, seq_str in strands:
            for frame in range(3):
                stops = []
                for i in range(frame, seq_len - 2, 3):
                    if seq_str[i:i+3] in Config.STOP_CODONS:
                        stops.append(i)
                
                if not stops:
                    continue
                    
                current_stop_idx = 0
                for i in range(frame, seq_len - 2, 3):
                    if seq_str[i:i+3] in Config.START_CODONS:
                        while current_stop_idx < len(stops) and stops[current_stop_idx] <= i:
                            current_stop_idx += 1
                            
                        if current_stop_idx < len(stops):
                            j = stops[current_stop_idx]
                            orf_len = j + 3 - i
                            
                            if orf_len >= min_orf_len:
                                prom_start = max(0, i - Config.MAX_PROMOTER_DIST)
                                prom_seq = seq_str[prom_start:i]
                                
                                reg_pssms = {k: v for k, v in pssms.items() if "RBS" not in k and "Sigma" not in k}
                                
                                # --- LA RÈGLE DES 25% ---
                                # On calcule tout à 25% minimum pour le mettre en cache
                                regs = scan_regulators(prom_seq, reg_pssms, 0.25)
                                rbs = find_rbs(prom_seq, pssms, 0.25)
                                proms = find_promoters(prom_seq, pssms, 0.25, _custom_spacers)
                                # ------------------------
                                
                                end_downstream = min(seq_len, j + 3 + Config.TERM_ZONE)
                                downstream_seq = seq_str[j+3 : end_downstream]
                                term = find_terminator(downstream_seq)
                                
                                orf_seq = seq_str[i:j+3]
                                gc3_val = calculate_gc3(orf_seq)
                                
                                try:
                                    prot_seq = str(Seq(orf_seq).translate(table=11, to_stop=True))
                                except Exception:
                                    prot_seq = "Translation Error"

                                pos_start = seq_len - i if strand_name == "Reverse (-)" else i + 1
                                pos_stop = seq_len - (j + 2) if strand_name == "Reverse (-)" else j + 3
                                    
                                all_raw_orfs.append({
                                    "Contig": contig_id, "Contig_Length": seq_len, "Strand": strand_name, 
                                    "Start": pos_start, "Stop": pos_stop, "Length": orf_len, 
                                    "Start_Codon": seq_str[i:i+3], "Regulators": regs, 
                                    "RBS": rbs, "Promoters": proms, "Terminator": term,
                                    "GC3": gc3_val, "Protein": prot_seq       
                                })
                                
    return all_raw_orfs, total_analyzed_len


def apply_dynamic_filters(raw_orfs, reg_thresh, prom_thresh, rbs_thresh):
    """ÉTAGE 2 : Le filtrage instantané via l'interface UI"""
    # 1. On travaille sur une copie pour ne pas détruire les données du cache
    working_orfs = copy.deepcopy(raw_orfs)
    all_final_operons = []
    
    # 2. On regroupe nos données brutes par contig
    orfs_by_contig = {}
    for orf in working_orfs:
        cid = orf['Contig']
        if cid not in orfs_by_contig:
            orfs_by_contig[cid] = []
        orfs_by_contig[cid].append(orf)
        
    for cid, contig_orfs in orfs_by_contig.items():
        
        # 3. ON APPLIQUE TES CURSEURS (Suppression instantanée de ce qui est sous le seuil)
        for orf in contig_orfs:
            orf['Regulators'] = [r for r in orf['Regulators'] if r['Percentage'] >= reg_thresh * 100]
            orf['Promoters'] = [p for p in orf['Promoters'] if p['Percentage'] >= prom_thresh * 100]
            orf['RBS'] = [r for r in orf['RBS'] if r['Percentage'] >= rbs_thresh * 100]

        # 4. On re-note les gènes avec ce qu'il leur reste
        max_found_len = max((orf['Length'] for orf in contig_orfs), default=1)
        for orf in contig_orfs:
            stars, numeric_score = assign_confidence_score(orf, max_found_len)
            orf['Stars'] = stars
            orf['Score_Cache'] = numeric_score
            
        # 5. On choisit les Champions (Filtrage des starts alternatifs)
        unique_orfs = {}
        for orf in contig_orfs:
            key = (orf['Strand'], orf['Stop'])
            if key not in unique_orfs:
                unique_orfs[key] = []
            unique_orfs[key].append(orf)
            
        found_orfs_clean = []
        for key, alternatives in unique_orfs.items():
            alternatives.sort(key=lambda x: (x['Score_Cache'], x['Length']), reverse=True)
            champion = alternatives[0]
            champion['Alternatives'] = alternatives[1:] 
            found_orfs_clean.append(champion)
            
        # 6. On regroupe les opérons
        contig_operons = group_operons(found_orfs_clean)
        
        for op in contig_operons:
            op['Contig'] = cid
            op['Contig_Length'] = contig_orfs[0]['Contig_Length']
            
        all_final_operons.extend(contig_operons)
        
    return all_final_operons

def group_operons(raw_orfs):
    orfs_plus = sorted([orf for orf in raw_orfs if orf['Strand'] == "Direct (+)"], key=lambda x: x['Start'])
    orfs_minus = sorted([orf for orf in raw_orfs if orf['Strand'] == "Reverse (-)"], key=lambda x: x['Start'], reverse=True)

    temp_groups = []
    for strand_orfs in [orfs_plus, orfs_minus]:
        if not strand_orfs: continue

        current_op = [strand_orfs[0]]
        for i in range(1, len(strand_orfs)):
            prev_gene = current_op[-1]
            curr_gene = strand_orfs[i]
            dist = curr_gene['Start'] - prev_gene['Stop'] if prev_gene['Strand'] == "Direct (+)" else prev_gene['Stop'] - curr_gene['Start'] 

            if -Config.MAX_OVERLAP <= dist <= Config.MAX_GAP:
                current_op.append(curr_gene)
            else:
                temp_groups.append(current_op)
                current_op = [curr_gene]

        if current_op:
            temp_groups.append(current_op)

    final_operons = []
    for group in temp_groups:
        avg_score = sum(g['Score_Cache'] for g in group) / len(group)
        if len(group) > 1 and len(group[0]['Promoters']) > 0 and group[-1]['Terminator'] is not None:
            avg_score += 10
            
        final_operons.append({"Genes": group, "Operon_Score": min(avg_score, 100)})

    # CORRECTION : Tri final par POSITION physique et non plus par note
    return sorted(final_operons, key=lambda x: min(g['Start'] for g in x['Genes']))

def extract_distance(dist_str):
    match = DIGIT_PATTERN.search(str(dist_str))
    return int(match.group()) if match else 0

def draw_operon_map(operon, total_len):
    features = []
    strand_val = 1 if operon[0]['Strand'] == "Direct (+)" else -1
    num_genes = len(operon)
    
    # --- LA CORRECTION DU ZOOM ---
    # On trouve les vraies limites de l'opéron
    toutes_les_positions = [pos for gene in operon for pos in (gene['Start'], gene['Stop'])]
    pos_min_op = min(toutes_les_positions)
    pos_max_op = max(toutes_les_positions)
    
    # On définit une fenêtre de vision avec 1000 pb de marge de chaque côté
    marge = 1000
    fenetre_debut = max(0, pos_min_op - marge)
    fenetre_fin = min(total_len, pos_max_op + marge)
    # ------------------------------
    
    for i, orf in enumerate(operon, 1):
        pos_min, pos_max = min(orf['Start'], orf['Stop']), max(orf['Start'], orf['Stop'])
        
        gene_color = "#4C72B0" if orf['Score_Cache'] > 70 else "#55A868" if orf['Score_Cache'] > 50 else "#C44E52"
        features.append(GraphicFeature(start=pos_min, end=pos_max, strand=strand_val, color=gene_color, label=f"Gene {i} ({orf['Length']} bp)"))
        
        if i == 1:
            if orf['Promoters']:
                prom = orf['Promoters'][0]
                dist = extract_distance(prom['Distance'])
                p_start = max(0, pos_min - dist) if strand_val == 1 else max(0, min(total_len, pos_max + dist) - 35)
                features.append(GraphicFeature(start=p_start, end=p_start + 35, strand=strand_val, color="#95a5a6", label=f"Prom {prom['Type']}"))

            for reg in orf['Regulators']:
                dist = extract_distance(reg['Position'])
                reg_len = len(reg['Sequence'])
                rg_start = max(0, pos_min - dist) if strand_val == 1 else max(0, min(total_len, pos_max + dist) - reg_len)
                features.append(GraphicFeature(start=rg_start, end=rg_start + reg_len, strand=strand_val, color="#9b59b6", label=reg['Regulator']))

        if orf['RBS']:
            dist = extract_distance(orf['RBS'][0]['Distance'])
            r_start = max(0, pos_min - dist) if strand_val == 1 else max(0, min(total_len, pos_max + dist) - 6)
            features.append(GraphicFeature(start=r_start, end=r_start + 6, strand=strand_val, color="#e67e22", label="RBS"))

        if i == num_genes:
            if orf['Terminator']:
                term = orf['Terminator']
                dist = extract_distance(term['Distance'])
                term_len = (len(term['Stem1']) * 2) + len(term['Loop']) + 6
                t_start = min(total_len, pos_max + dist) if strand_val == 1 else max(0, pos_min - dist - term_len)
                features.append(GraphicFeature(start=t_start, end=t_start + term_len, strand=strand_val, color="#c0392b", label="Term"))

    # On force le dessin à ne prendre en compte que notre fenêtre zoomée
    record = GraphicRecord(first_index=fenetre_debut, sequence_length=fenetre_fin, features=features)
    fig, ax = plt.subplots(figsize=(11, 3.5))
    record.plot(ax=ax)
    
    # On cadre l'axe X strictement sur notre fenêtre
    ax.set_xlim(fenetre_debut, fenetre_fin)
    
    # On recalcule les repères (ticks) pour qu'ils soient lisibles sur cette petite échelle
    taille_fenetre = fenetre_fin - fenetre_debut
    pas_tick = max(100, taille_fenetre // 10) 
    ax.set_xticks(list(range(fenetre_debut, fenetre_fin + 1, pas_tick)))
    
    ax.tick_params(axis='x', labelsize=8)
    fig.patch.set_alpha(0.0) 
    return fig

# ==========================================
# 4. EXPORT & UTILS
# ==========================================

def generate_excel_export(operon_results):
    rows = []
    operons_by_contig = {}
    
    for op in operon_results:
        cid = op.get('Contig', 'Unknown')
        if cid not in operons_by_contig:
            operons_by_contig[cid] = []
        operons_by_contig[cid].append(op)
        
    for cid, op_list in operons_by_contig.items():
        for index_op, op in enumerate(op_list, 1):
            for i_gene, orf in enumerate(op['Genes'], 1):
                
                prom_info = f"[{orf['Promoters'][0]['Type']}] Aff: {orf['Promoters'][0]['Percentage']}%" if orf['Promoters'] else "None"
                rbs_info = f"Seq: {orf['RBS'][0]['Sequence']} | Aff: {orf['RBS'][0]['Percentage']}%" if orf['RBS'] else "None"
                reg_info = " / ".join([f"{r['Regulator']} ({r['Percentage']}%)" for r in orf['Regulators']]) if orf['Regulators'] else "None"
                term_info = f"Thermo: {orf['Terminator']['Thermo_Score']}" if orf['Terminator'] else "None"

                pos_min, pos_max = min(orf['Start'], orf['Stop']), max(orf['Start'], orf['Stop'])

                rows.append({
                    "Contig": cid, "Operon ID": f"Operon {index_op}", "Operon Score (/100)": round(op['Operon_Score'], 1),
                    "Gene ID": f"Gene {i_gene}", "Strand": orf['Strand'], "Start": pos_min, "Stop": pos_max,
                    "Length (bp)": orf['Length'], "Gene Score (/100)": round(orf['Score_Cache'], 1),
                    "GC3 (%)": round(orf['GC3'], 1), "Main Promoter": prom_info,
                    "Main RBS": rbs_info, "Regulators": reg_info, "Terminator": term_info,
                    "Protein Sequence": orf['Protein']
                })
                
    buffer = BytesIO()
    with pd.ExcelWriter(buffer, engine='openpyxl') as writer:
        pd.DataFrame(rows).to_excel(writer, index=False, sheet_name='M.A.P.P.E.R Results')
    return buffer.getvalue()

def generate_template_bytes():
    prom_data, reg_data = [], []
    for sig in ["SigmaA", "SigmaB", "SigmaE", "SigmaK", "SigmaG"]:
        for s35, s10 in zip_longest(DEFAULT_REGULATORS.get(f"{sig}_-35", []), DEFAULT_REGULATORS.get(f"{sig}_-10", []), fillvalue=""):
            prom_data.append({"Factor_Name": sig, "Box_-35": s35, "Box_-10": s10, "Spacer": "13-15" if sig in ["SigmaB", "SigmaE", "SigmaK"] else "15-18"})

    for reg in ["CodY", "PlcR", "CRE high affinity", "CRE low affinity", "CRE activator", "Spo0A"]:
        for seq in DEFAULT_REGULATORS.get(reg, []):
            reg_data.append({"Regulator_Name": reg, "Sequence": seq})

    rbs_data = [{"Sequence": s} for s in DEFAULT_REGULATORS.get("RBS (Shine-Dalgarno)", [])]

    buffer = BytesIO()
    with pd.ExcelWriter(buffer, engine='openpyxl') as writer:
        pd.DataFrame(prom_data).to_excel(writer, sheet_name='Promoters', index=False)
        pd.DataFrame(reg_data).to_excel(writer, sheet_name='Regulators', index=False)
        pd.DataFrame(rbs_data).to_excel(writer, sheet_name='RBS', index=False)
    return buffer.getvalue()

# ==========================================
# 5. STREAMLIT USER INTERFACE
# ==========================================

def main_ui():
    st.set_page_config(page_title="M.A.P.P.E.R Sequence Analyzer", page_icon="🧬", layout="wide")
    st.title("🧬 M.A.P.P.E.R")
    st.markdown("ORF, Promoter, and Regulator Search")
    st.divider()

    with st.sidebar:
        st.header("⚙️ Analysis Parameters")
        # CORRECTION MAJEURE ICI : step=1
        min_orf_len = st.number_input("Minimum gene size (bp)", min_value=50, max_value=5000, value=500, step=1)
        
        st.subheader("Detection Thresholds")
        reg_thresh = st.slider("Regulator Affinity", min_value=0.0, max_value=1.0, value=0.70, step=0.05)
        prom_thresh = st.slider("Promoter Affinity", min_value=0.0, max_value=1.0, value=0.70, step=0.05)
        rbs_thresh = st.slider("RBS Affinity", min_value=0.0, max_value=1.0, value=0.65, step=0.05)
        
        st.divider()
        st.subheader("🧫 Affinity Matrices")
        
        st.download_button(
            label="📄 Download Excel Template",
            data=generate_template_bytes(),
            file_name="Template_Matrices_MAPPER.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            help="Download this file, edit for your organism, and upload below!"
        )
        custom_matrices_file = st.file_uploader("Optional: Load a dictionary (.xlsx)", type=["xlsx"])

    uploaded_file = st.file_uploader("📥 Upload your DNA sequence (.fasta, .fa, .fna)", type=["fasta", "fa", "fna"])
    
    if uploaded_file:
        if custom_matrices_file:
            active_dict, active_spacers = read_custom_matrices(custom_matrices_file)
            st.info("💡 Using custom multi-sheet dictionary.")
        else:
            active_dict, active_spacers = DEFAULT_REGULATORS, DEFAULT_SPACERS
            st.info("💡 Using default dictionary (Bacillus).")
            
        seq_text = uploaded_file.getvalue().decode("utf-8")
        
        # TEMPS 1 : L'extraction lourde mise en cache (Ne dépend plus des curseurs)
        with st.spinner("🧬 Extracting Genomic Data (Runs once)..."):
            raw_results_tuple = extract_raw_orfs(
                seq_text, min_orf_len, active_dict, active_spacers
            )
            
        if raw_results_tuple:
            raw_orfs, total_len = raw_results_tuple
            
            # TEMPS 2 : Le filtrage instantané (Réagit aux curseurs)
            with st.spinner("⚡ Applying Dynamic Filters..."):
                orf_results = apply_dynamic_filters(raw_orfs, reg_thresh, prom_thresh, rbs_thresh)
                st.success(f"Dynamic analysis complete! {len(orf_results)} operons/genes found across {total_len} bp.")
            
            st.download_button(
                label="📥 Download full table (Excel)",
                data=generate_excel_export(orf_results),
                file_name="Results_MAPPER.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                type="primary" 
            )
            st.divider()
            
            operons_by_contig = {}
            for op in orf_results:
                cid = op.get('Contig', 'Unknown')
                if cid not in operons_by_contig: operons_by_contig[cid] = []
                operons_by_contig[cid].append(op)
            
            if not operons_by_contig:
                st.warning("⚠️ No genes found matching these parameters.")
            else:
                list_contigs = list(operons_by_contig.keys())
                
                if len(list_contigs) > 1:
                    st.markdown("### 📂 Multi-Fragment Navigation")
                    selected_contig = st.selectbox("Select the fragment to display:", list_contigs)
                    # On ne garde que le contig sélectionné pour l'affichage
                    contigs_to_display = {selected_contig: operons_by_contig[selected_contig]}
                else:
                    contigs_to_display = operons_by_contig

                for contig_id, op_list in contigs_to_display.items():
                    c_len = op_list[0].get('Contig_Length', total_len)
                    
                    st.write("") 
                    st.markdown(f"## 🧩 Fragment: `{contig_id}`")
                    st.caption(f"📏 Size: **{c_len} bp** | 🔗 Found Blocks: **{len(op_list)}** | 🧬 Total ORFs: **{sum(len(o['Genes']) for o in op_list)}**")
                    
                    # --- LE RESTE DE TON CODE NE CHANGE PAS (Filtres, Pagination, etc.) ---
                    col_sort, col_page = st.columns([2, 1])
                    with col_sort:
                        sort_choice = st.selectbox(
                            "Sort results by:",
                            [
                                "📍 Position (Left to Right)", 
                                "⭐ Operon Score (Descending)", 
                                "📏 Total Size (Descending)", 
                                "🔗 Number of Genes (Descending)"
                            ],
                            key=f"sort_{contig_id}"
                        )
                    
                # Application logique du tri choisi
                if "Position" in sort_choice:
                    op_list.sort(key=lambda op: min(g['Start'] for g in op['Genes']))
                elif "Score" in sort_choice:
                    op_list.sort(key=lambda op: op['Operon_Score'], reverse=True)
                elif "Size" in sort_choice:
                    op_list.sort(key=lambda op: max(g['Stop'] for g in op['Genes']) - min(g['Start'] for g in op['Genes']), reverse=True)
                elif "Number" in sort_choice:
                    op_list.sort(key=lambda op: len(op['Genes']), reverse=True)

                # --- NOUVEAUTÉ 2 : LA PAGINATION (LE SAUVEUR DE RAM) ---
                ITEMS_PER_PAGE = 15
                total_pages = max(1, (len(op_list) - 1) // ITEMS_PER_PAGE + 1)
                
                with col_page:
                    if total_pages > 1:
                        current_page = st.number_input(f"Page (1 - {total_pages})", min_value=1, max_value=total_pages, value=1, step=1, key=f"page_{contig_id}")
                    else:
                        current_page = 1
                        st.info("All results displayed.")
                        
                start_idx = (current_page - 1) * ITEMS_PER_PAGE
                end_idx = start_idx + ITEMS_PER_PAGE
                
                st.divider()
                
                # On boucle UNIQUEMENT sur les 15 éléments de la page actuelle !
                for index_op, data_op in enumerate(op_list[start_idx:end_idx], start=start_idx + 1):
                    operon = data_op['Genes']
                    brin_op = operon[0]['Strand']
                    toutes_les_positions = [pos for gene in operon for pos in (gene['Start'], gene['Stop'])]
                    op_len = max(toutes_les_positions) - min(toutes_les_positions) + 1
                    
                    titre = f"🧬 ISOLATED GENE [{index_op}]" if len(operon) == 1 else f"🔗 OPERON [{index_op}] | {len(operon)} genes"
                    titre += f" | Strand {brin_op} | Score: {data_op['Operon_Score']:.1f}/100 | {op_len} bp"
                    
                    with st.expander(titre):
                        fig_carte = draw_operon_map(operon, c_len)
                        st.pyplot(fig_carte)
                        plt.close(fig_carte)
                        st.divider() 
                        
                        for i_gene, orf in enumerate(operon, 1):
                            if len(operon) > 1: st.markdown(f"### 🔹 Gene {i_gene} (Score: {orf['Score_Cache']:.1f}/100)")
                            
                            st.markdown(f"**Position:** {orf['Start']} ➡️ {orf['Stop']} (Start: `{orf['Start_Codon']}`) | **Length:** `{orf['Length']} bp`")
                            
                            if orf.get('Alternatives'):
                                with st.expander("🔄 Alternative Start Codons (Primer Check)"):
                                    for alt in orf['Alternatives']:
                                        st.write(f"- Start at **{alt['Start']}** (`{alt['Start_Codon']}`) | Length: {alt['Length']} bp | Score: {alt['Score_Cache']:.1f}/100")
                            
                            col1, col2 = st.columns(2)
                            with col1:
                                if i_gene == 1:
                                    st.markdown("#### 🟢 Promoters")
                                    if orf['Promoters']:
                                        for p in orf['Promoters']: st.code(f"[{p['Type']}] Pos: {p['Distance']:<7} | Aff: {p['Percentage']:.1f}%\n[-35] {p['Seq_35']} --({p['Spacer']}bp)-- [-10] {p['Seq_10']}")
                                    else: st.write("None")
                                st.markdown("#### 🔵 RBS")
                                if orf['RBS']:
                                    for r in orf['RBS']: st.code(f"Seq: {r['Sequence']} | Pos: {r['Distance']:<7} | Aff: {r['Percentage']:.1f}%")
                                else: st.write("None")
                                    
                            with col2:
                                if i_gene == 1:
                                    st.markdown("#### 🟠 Regulators")
                                    if orf['Regulators']:
                                        for r in sorted(orf['Regulators'], key=lambda x: x['Score'], reverse=True):
                                            st.code(f"{r['Regulator']} | Aff: {r['Percentage']:.1f}%\nSeq: {r['Sequence']} | Pos: {r['Position']}")
                                    else: st.write("None")
                                if i_gene == len(operon):
                                    st.markdown("#### 🔴 Terminator")
                                    if orf['Terminator']:
                                        t = orf['Terminator']
                                        st.code(f"Pos: {t['Distance']} | Thermo: {t['Thermo_Score']}\n[{t['Stem1']}]-{t['Loop']}-[{t['Stem2']}]-{t['Tail']}...")
                                    else: st.write("None")
                            
                            with st.expander("🧬 View translated protein sequence"): st.code(orf['Protein'])
                            if len(operon) > 1 and i_gene < len(operon): st.divider()

if __name__ == "__main__":
    main_ui()