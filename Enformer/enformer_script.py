import re, csv, kipoiseq, requests, tensorflow_hub as hub, tensorflow as tf
from pathlib import Path
from tqdm import tqdm

print("Loading Enformer …")
model = hub.load("https://tfhub.dev/deepmind/enformer/1").model
print("Model ready.\n")

FASTA_IN      = "/project2/kribelba_1515/saadawy/Enformer/upstream/enformer_inputs_upstream.fasta"
OUT_TSV       = "enformer_upstream/enformer_predictions_unfiltered.tsv"

LEFT_HOM      = "ATAACTTCGTATAGGATACTTTATACGAAGTTAT"
RIGHT_HOM     = "TAACTTCGTATAATGTATGCTATACGAAGTTAT"
ENH_PREFIX    = "CGTGAGAGAACGCTC"        
ENH_SUFFIX    = "ACGGATCCGACAGC"      
ENH_OFFSET    = 40                     # enhancer start within construct
SEQ_LEN       = 393_216
ENFORMER_PRED_WINDOW = 114_688
MID_START     = SEQ_LEN//2 - ENFORMER_PRED_WINDOW//2
BIN_SIZE      = 128
TRACKS        = [5110, 69, 688]         # CAGE, DNase, H3K27ac

def fasta_iter(fp):
    with open(fp) as fh:
        hdr, seq = None, []
        for ln in fh:
            if ln.startswith(">"):
                if hdr:
                    yield hdr, "".join(seq)
                hdr, seq = ln[1:].rstrip(), []
            else:
                seq.append(ln.rstrip())
        if hdr:
            yield hdr, "".join(seq)

def classify_bin(bin_start, bin_end,
                 c_start, c_end,
                 enh_start, enh_end):
    """Return a string describing what the 128-bp bin overlaps."""
    tags = []

    # whole bin left of construct
    if bin_end   < c_start:  return "left_flank"
    # whole bin right of construct
    if bin_start > c_end:    return "right_flank"

    # overlaps construct
    if bin_end   < enh_start:
        tags.append("enhancer_upstream")
    elif bin_start > enh_end:
        tags.append("enhancer_downstream")
    else:
        tags.append("enhancer")

    # might also overlap flank if bin straddles boundary
    if bin_start < c_start:
        tags.append("left_flank")
    if bin_end   > c_end:
        tags.append("right_flank")

    return " + ".join(tags)
    
total_seqs = sum(1 for ln in open(FASTA_IN) if ln.startswith(">"))
print(f"Processing {total_seqs} sequences …\n")

with open(OUT_TSV, "w", newline="") as out_fh:
    writer = csv.writer(out_fh, delimiter="\t")
    writer.writerow(["Enhancer_ID", "Bin", "Shift", "Bin_Info", "CAGE", "DNase", "H3K27ac"])

    for hdr, seq in tqdm(fasta_iter(FASTA_IN),
                         total=total_seqs,
                         ncols=80,
                         desc="Enformer predictions"):

        # ---------- parse header -----------------------------------
        enh_id, rest = hdr.split("_posShift")
        shift_val    = int(rest.split("_")[0])          # -58, 6, …

        # ---------- locate construct & enhancer --------------------
        c_start = seq.find(LEFT_HOM)
        c_end   = seq.find(RIGHT_HOM) + len(RIGHT_HOM)

        enh_start = seq.find(ENH_PREFIX, c_start)
        enh_end   = seq.find(ENH_SUFFIX, enh_start) + len(ENH_SUFFIX) - 1

        # ---------- Enformer prediction ----------------------------
        one_hot = kipoiseq.transforms.functional.one_hot_dna(seq)\
                                                    .astype("float32")[None, ...]
        preds = model.predict_on_batch(one_hot)['human'][0].numpy()[:, TRACKS]

        # ---------- write 896 rows ---------------------------------
        for b in range(896):
            bin_start = MID_START + b*BIN_SIZE
            bin_end   = bin_start + BIN_SIZE - 1
            info = classify_bin(bin_start, bin_end,
                                c_start, c_end,
                                enh_start, enh_end)

            writer.writerow([enh_id, f"{b:03}", shift_val, info,
                             preds[b,0], preds[b,1], preds[b,2]])

print("✅  Finished TSV:", OUT_TSV)
