import pandas as pd
import argparse

def read_gtf(file):
    gtf_data = []
    with open(file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            start = int(parts[3])
            end = int(parts[4])
            attributes = parts[8]
            gtf_data.append([parts[0], parts[1], parts[2], start, end, parts[5], parts[6], parts[7], attributes])
    return pd.DataFrame(gtf_data, columns=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes'])

def overlap_percentage(start1, end1, start2, end2):
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    overlap = max(0, overlap_end - overlap_start + 1)
    len1 = end1 - start1 + 1
    len2 = end2 - start2 + 1
    return overlap / len1, overlap / len2

def filter_overlapping_annotations(gtf_df):
    to_remove = set()
    gtf_df = gtf_df.sort_values(by=['seqname', 'start', 'end']).reset_index(drop=True)

    for i in range(len(gtf_df) - 1):
        for j in range(i + 1, len(gtf_df)):
            if gtf_df.at[i, 'seqname'] != gtf_df.at[j, 'seqname']:
                continue
            if gtf_df.at[j, 'start'] > gtf_df.at[i, 'end']:
                break
            overlap1, overlap2 = overlap_percentage(gtf_df.at[i, 'start'], gtf_df.at[i, 'end'],
                                                    gtf_df.at[j, 'start'], gtf_df.at[j, 'end'])
            if overlap1 > 0.1 or overlap2 > 0.1:
                if (gtf_df.at[i, 'end'] - gtf_df.at[i, 'start']) <= (gtf_df.at[j, 'end'] - gtf_df.at[j, 'start']):
                    to_remove.add(i)
                else:
                    to_remove.add(j)

    filtered_gtf_df = gtf_df.drop(list(to_remove)).reset_index(drop=True)
    return filtered_gtf_df

def write_gtf(df, file):
    with open(file, 'w') as f:
        for index, row in df.iterrows():
            f.write(f"{row['seqname']}\t{row['source']}\t{row['feature']}\t{row['start']}\t{row['end']}\t{row['score']}\t{row['strand']}\t{row['frame']}\t{row['attributes']}\n")

if __name__ == "__main__":
    # Set up argument parsing
    parser = argparse.ArgumentParser(description='Filter overlapping GTF annotations.')
    parser.add_argument('input_file', type=str, help='Input GTF file')
    parser.add_argument('output_file', type=str, help='Output GTF file')
    
    # Parse the arguments
    args = parser.parse_args()

    # Process the GTF file
    gtf_df = read_gtf(args.input_file)
    filtered_gtf_df = filter_overlapping_annotations(gtf_df)
    write_gtf(filtered_gtf_df, args.output_file)
