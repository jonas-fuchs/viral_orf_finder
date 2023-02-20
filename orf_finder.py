from Bio import SeqIO
import pandas as pd
import warnings
import argparse
import sys

def extract_ref_info(path_to_gb, to_extract = ["CDS", "source"]):
    """
    extract user defined information from genbank file as dictionary
    """
    gb_file = SeqIO.read(path_to_gb, format="genbank")
    # initialize stuff
    fasta = gb_file.seq
    feature_list = []
    final_dict = dict()
    # extract relevant infos from *.gb file
    for ann in to_extract:
        annotation_dict = dict()
        counter = 0 #countes how often NA has to be written
        for features in gb_file.features:
            feature_list.append(features.type)
            if features.type == ann:
                items_temp = dict(features.qualifiers.items())
                items_temp["location"] = [min(features.location),max(features.location)+1]
                items_temp["seq"] = [str(fasta[items_temp["location"][0]:items_temp["location"][1]])]
                for entry in items_temp:
                    if entry in annotation_dict:
                        while len(annotation_dict[entry]) < counter:
                            annotation_dict[entry].append("NA")
                    else:
                        annotation_dict[entry] = counter*["NA"]
                    if len(items_temp[entry]) > 1:
                        annotation_dict[entry].append(items_temp[entry])
                    else:
                        annotation_dict[entry].append(items_temp[entry][0])
                counter += 1
        if ann not in feature_list:
            print(ann, " was not found in the provided *.gb file")
        else:
            final_dict[ann]=annotation_dict
    return final_dict

def get_params(path_to_gb):
    """
    extract parameters for orf finder
    """
    location_CDS = extract_ref_info(path_to_gb, to_extract = ["CDS"])["CDS"]["location"]
    previous_start = 0
    previous_stop = 0
    length = []
    no_overlap = []
    for location in location_CDS:
        length.append(location[1]-location[0])
        if location[1] > previous_stop:
            no_overlap.append(True)
        elif location[0] == previous_start:
            no_overlap.append(True)
        else:
            no_overlap.append(False)
        previous_start = location[0]
        previous_stop = location[1]
    if any(item is False for item in no_overlap):
        no_overlap = False
    else:
        no_overlap = True
    params = (min(length), no_overlap)

    return params

def find_orfs(path_to_fasta, min_len, start_codons = ["ATG"], stop_codons = ["TAG", "TGA", "TAA"], internal = False, circular = False, partial = False, no_overlap = False, strands = "+,-"):

    # translation function results in an error if seq is not a multiple of three
    warnings.filterwarnings('ignore')

    """
    - finds ORFs and their translations -

    possible ORF types:
    0: -----[M-------------*]-------- complete
    1: ------M--[M---------*]-------- complete_internal
    2: -----[M----------------------- 5_partial
    3: ------M--[M------------------- 5_partial_internal
    4: [--------------------*]------- 3_partial
    5: [----------------------------] 5_3_partial
    6: ------*]-------------[M------- circular
    7: ------*]--------------M---[M-- circular_internal

    no overlap algorithm:
    frame 1: -[M------*]--- ----[M--*]---------[M-----
    frame 2: -------[M------*]---------[M---*]--------
    frame 3: [M---*]-----[M----------*]----------[M---

    results: [M---*][M------*]--[M--*]-[M---*]-[M-----
    frame:    3      2           1      2       1
    """
    # for no overlap internal orfs have to be considered
    if no_overlap:
        internal = True
    # initialize
    final_orfs = []
    orfs_written = 0
    # read in fasta
    fasta = SeqIO.read(path_to_fasta, format="fasta").seq
    fasta_name = SeqIO.read(path_to_fasta, format="fasta").id
    # sanity check: convert to upper and transcribe
    fasta = fasta.upper()
    if "U" in fasta:
        fasta = fasta.back_transcribe()

    # define which strands to look at
    separators = ["; ",", "," ","\t",";",","," ","\t"]
    for sep in separators:
        if sep in strands:
            strands = strands.split(sep)
            strands = set(strands)
            break
    else:
        strands = set([strands])

    if "+" in strands and "-" in strands:
        strand_list = [("+", fasta), ("-", fasta.reverse_complement())]
    elif "+" in strands:
        strand_list = [("+", fasta)]
    elif "-" in strands:
        strand_list = [("-", fasta.reverse_complement())]
    else:
        sys.exit("strand must be -,+")
    # iterate over strand and complement
    for strand, seq in strand_list:
        seq = str(seq)
        if circular:
            seq_circ = seq+seq
        for frame in range(3):
            # track orfs with [start, stop, ORF_type]
            orf = [-1,-1,-1]
            internal_starts = []
            orf_list = []
            starts_and_stops = []
            count_stops = 0
            # for each frame remember [type, position]
            # type:  0 = start, 1 = stop, 2 = dummy_start, 3 = dummy stop
            for i in range(frame, len(seq),3):
                found_stop = "no"
                # search for start codons
                if i < len(seq)-2:
                    if any(seq[i:i+3] in s for s in start_codons):
                        starts_and_stops.append([0,i])
                    # create dummy start if start of the seq
                    # is no start codon
                    elif i == frame:
                        starts_and_stops.append([2,frame])
                    # search for stop codons
                    if any (seq[i:i+3] in s for s in stop_codons):
                        starts_and_stops.append([1,i+3])
                        count_stops += 1
                        found_stop = "yes"
                # create dummy stop
                if i >= len(seq)-3 and found_stop == "no":
                    starts_and_stops.append([3,len(seq)])

            # search for orfs
            # if the codon after the dummy start is a stop
            # edge case: 3_partial
            if partial:
                if starts_and_stops[0][0] == 2 and starts_and_stops[1][0] == 1:
                    orf = [frame, starts_and_stops[1][1],4]
                    if orf[1]-orf[0] >= min_len:
                        orf_list.append(orf)
                    # start new orf
                    orf = [-1,-1,-1]
            # if no stops have been counted
            # edge case: 5_3_partial
            if count_stops == 0:
                if partial:
                    orf = [frame, len(seq),5]
                    orf_list.append(orf)
                    # start new orf
                    orf = [-1,-1,-1]
            # search for all other orfs
            else:
                for position in starts_and_stops:
                    # no start has been written so far and position is a start
                    # -> start of an orf
                    if orf[0] == -1 and position[0] == 0:
                        orf[0] = position[1]
                    # start has been found, its not the already written start, no stop
                    # was written and position is a start
                    # -> internal orf
                    if internal:
                        if all((orf[0] != -1, orf[0] != position[1], orf[1] == -1, position[0] == 0)):
                            # remember all positions of internal starts
                            internal_starts.append(position[1])
                    # start has been found and position is a stop
                    # -> complete orf or orf continues outside sequence (dummy stop)
                    if orf[0] != -1 and any((position[0] == 1, position[0] == 3)):
                        orf[1] = position[1]
                        # real stop
                        if position[0] == 1:
                            orf[2] = 0
                        # edge case - reached dummy stop, orf continues or is circular
                        elif position[0] == 3:
                            # orf continues
                            if partial:
                                orf[2] = 2
                            else:
                                orf = [-1,-1,-1]
                            # consider a circular orf
                            if circular:
                                orf_circ = [orf[0], -1, 6]
                                # iterate again over apended sequence
                                # and stop at first stop codon
                                for i in range(orf_circ[0], len(seq_circ),3):
                                    if any(seq_circ[i:i+3] in s for s in stop_codons):
                                        orf_circ[1] = i+3 - len(seq)
                                        break
                                if len(seq)-orf_circ[0]+orf_circ[1] >= min_len and orf[1] != -1:
                                    orf_list.append(orf_circ)
                                    # consider internal circular orfs
                                    if internal:
                                        for start in internal_starts:
                                            if len(seq)-start+orf_circ[1] >= min_len:
                                                internal_orf_circ = [start, orf_circ[1], 7]
                                                orf_list.append(internal_orf_circ)
                        if orf[1]-orf[0] >= min_len:
                            orf_list.append(orf)
                            # a orf with the right len was found
                            # write all internal orfs
                            if internal:
                                for start in internal_starts:
                                    if orf[1]-start >= min_len:
                                        if position[0] == 1:
                                            internal_orf = [start, orf[1], 1]
                                        # edge case - continues outside and internal
                                        elif position[0] == 3:
                                            internal_orf = [start, orf[1], 3]
                                        orf_list.append(internal_orf)
                        # start new orf
                        orf = [-1,-1,-1]
                        # start new list of internal starts
                        internal_starts = []

            for orf in orf_list:
                # edge case: circular and circular_internal
                if any((orf[2] == 6, orf[2] == 7)):
                    fasta_to_trans = fasta+fasta
                else:
                    fasta_to_trans = fasta
                if strand == "-":
                    # get the right coordniates on the original sequence
                    orf[0] = len(seq) - orf[0]
                    orf[1] = len(seq) - orf[1]
                    orf = [orf[1],orf[0], orf[2]]
                    # edge case: circular and circular_internal
                    if any((orf[2] == 6, orf[2] == 7)):
                        trans = fasta_to_trans[orf[0]:len(seq)+orf[1]].reverse_complement().translate()
                    else:
                        trans = fasta_to_trans[orf[0]:orf[1]].reverse_complement().translate()
                else:
                    # edge case: circular and circular_internal
                    if any((orf[2] == 6, orf[2] == 7)):
                        trans = fasta_to_trans[orf[0]:len(seq)+orf[1]].translate()
                    else:
                        trans = fasta_to_trans[orf[0]:orf[1]].translate()
                if orf[2] == 0:
                    status = "complete"
                elif orf[2] == 1:
                    status = "internal"
                elif orf[2] == 2:
                    status = "5_partial"
                elif orf[2] == 3:
                    status = "5_partial_internal"
                elif orf[2] == 4:
                    status = "3_partial"
                elif orf[2] == 5:
                    status = "5_3_partial"
                elif orf[2] == 6:
                    status = "circular"
                elif orf[2] == 7:
                    status = "circular_internal"
                final_orfs.append((fasta_name,"ORF_"+str(orfs_written),status, strand,
                frame+1, orf[0], orf[1], len(trans), str(trans)))
                orfs_written += 1
    # create dataframe
    final_orfs = pd.DataFrame(final_orfs, columns=["seq_name","orf_name", "type", "strand",
    "frame", "start","stop","as_length", "translation"])


    # creedy, checks for overlaps:
    # walks along the sequence until first start is found and ignores all other
    # starts until the stop is reached. checks also for the reverse strand.
    # for circular orfs starts are only ignored until the end of the sequence.
    if no_overlap:
        rows_to_keep = []
        for strand in strands:
            subset_df = final_orfs[final_orfs["strand"] == strand]
            if strand == "-":
                subset_df["start"] = len(seq) - subset_df["start"]
                subset_df["stop"] = len(seq) - subset_df["stop"]
                subset_df[["start", "stop"]] = subset_df[["stop", "start"]]
            subset_df = subset_df.sort_values(["start"])
            # if first orf is circular set stop to len(seq)
            if subset_df.iloc[0]["type"] == "circular":
                last_stop = len(seq)
            # otherwise remember the stop of the first orf
            else:
                last_stop = subset_df.iloc[0]["stop"]
            rows_to_keep.append(subset_df.index.values[0])
            # check if following starts are > last orf stop
            for i in range(1, len(subset_df)):
                if last_stop < subset_df.iloc[i]["start"]:
                    rows_to_keep.append(subset_df.index.values[i])
                    if subset_df.iloc[i]["type"] == "circular":
                        last_stop = len(seq)
                    else:
                        last_stop = subset_df.iloc[i]["stop"]
        final_orfs = final_orfs.iloc[rows_to_keep]
        # internal with no overlaps will get flagged as complete
        final_orfs = final_orfs.replace("internal", "complete")

    # sort by strand and starts and rename
    final_orfs = final_orfs.sort_values(["strand", "start"])
    final_orfs["orf_name"] = ["orf_"+str(n) for n in range(0, len(final_orfs))]

    return final_orfs




if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("infile", help="Filepath to fasta")
    parser.add_argument("-r", "--reference", type = str)
    parser.add_argument("-m", "--min_len", type = int, default=200, help="internal orfs")
    parser.add_argument("-i", "--internal", type = bool, default=True, help="internal orfs")
    parser.add_argument("-c", "--circular", type = bool, default=False, help="circular orfs")
    parser.add_argument("-p", "--partial", type = bool, default=False, help="partial orfs missing start and/or stop")
    parser.add_argument("-n", "--no-overlap", type = bool, default=False, help="search orfs with no overlap")
    parser.add_argument("-s", "--strands", type = str, default="+,-", help="+,- strand")

    args = parser.parse_args()


    if args.reference is not None:
        params = get_params(args.reference)
        print(
            find_orfs(
                args.infile,
                params[0],
                internal=args.internal,
                circular=args.circular,
                partial=args.partial,
                no_overlap = params[1],
                strands= args.strands
                )
            )
    else:
        print(
            find_orfs(
                args.infile,
                args.min_len,
                internal=args.internal,
                circular=args.circular,
                partial=args.partial,
                no_overlap = args.no_overlap,
                strands= args.strands
                )
            )
