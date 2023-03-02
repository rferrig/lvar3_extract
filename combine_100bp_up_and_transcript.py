from Bio import SeqIO
import re
from Bio.SeqRecord import SeqRecord
import pandas as pd
import pathlib

id_to_name = pd.read_csv("/projectnb/bradham/workflows/lvar3-0_extraction/top_probe_models.csv")
id_to_name.seq_id = id_to_name.seq_id.str.replace("(-RA$)", "", regex=True)
id_to_name = id_to_name.set_index("seq_id").to_dict()["name"]
offset_regex = re.compile("(?<=offset\:)[0-9]+(?!' ')")
transcripts = {
    x.id.split("::")[0].replace("-RA", ""): x
    for x in SeqIO.parse("regulons_transcripts.fa", "fasta")
}
upstream = {x.id.split("::")[0]: x for x in SeqIO.parse("100bp_up.fa", "fasta")}
downstream = {x.id.split("::")[0]: x for x in SeqIO.parse("100bp_down.fa", "fasta")}


def find_offset(x):
    return int(offset_regex.search(x).group())


def mark_stop(seq):
    stops = ["TAA", "TGA", "TAG"]
    i = 0
    found = False
    while i < len(seq) or found == False:
        if seq[i : i + 3] in stops:
            found = True
            break
        i += 3
    i = min(i, len(seq))
    return seq[:i] + seq[i:].lower()


outdir = pathlib.Path("output")
new_sequences = [None] * len(transcripts)
#print(transcripts.keys())
for i, seq_id in enumerate(transcripts.keys()):
    #print(seq_id)
    print(id_to_name[seq_id])
    offset = find_offset(transcripts[seq_id].description)
    if offset < 100:
        up_seq = upstream[seq_id][-(100 - offset) :].seq.lower()
    else:
        up_seq = transcripts[seq_id][(offset - 100) : offset].seq.lower()
    down_seq = transcripts[seq_id].seq.lower()
    new_seq = (
        up_seq
        + mark_stop(transcripts[seq_id].seq[offset:])
        + downstream[seq_id].seq.lower()
    )

    new_sequences[i] = SeqRecord(
        new_seq,
        id=upstream[seq_id].id.replace(seq_id, id_to_name[seq_id]),
        name=upstream[seq_id].id.replace(seq_id, id_to_name[seq_id]),
        description=f": {transcripts[seq_id].description.split(' ')[0]}".replace(
            seq_id, id_to_name[seq_id]
        ),
    )
    SeqIO.write(new_sequences[i], outdir.joinpath(id_to_name[seq_id]).with_suffix(".fa"), "fasta")

SeqIO.write(new_sequences, "regulons_100bp-up_TSS.fasta", "fasta")
