# Example Python script: extract_gsm_sra_from_hdf5.py
import h5py
import sys

hdf5_file_path = '/scratch/dhan/order_of_splicing/data/human_gene_perturbation_GSE_candidates_sample_metadata_hek293.hdf5'
output_tsv_path = 'gsm_sra_list_for_pipeline_hek293.tsv' # Output file

with h5py.File(hdf5_file_path, 'r') as hf, open(output_tsv_path, 'w') as outfile:
    if 'GSE_sample_info_strings' in hf:
        dataset = hf['GSE_sample_info_strings']
        print(f"Processing {dataset.shape[0]} entries from HDF5.")
        outfile.write("gsm_id\tsra_ids\n") # Header for the new TSV

        for i in range(dataset.shape[0]):
            entry_bytes = dataset[i]
            try:
                entry_str = entry_bytes.decode('utf-8')
            except AttributeError: # Already a string
                entry_str = entry_bytes

            lines = entry_str.strip().split('\n')
            # Series info is in lines[0] to lines[3]
            # Sample info starts from lines[4]
            for sample_line in lines[4:]:
                parts = sample_line.split('\t')
                if len(parts) == 3:
                    gsm_id = parts[0]
                    sra_ids = parts[2] # This should be your comma-separated SRRs or single SRX

                    # If SRA ID is an SRX, you might need to resolve it to SRRs here
                    # For now, we assume it's SRRs or the bash script handles SRX (less ideal)
                    # The bash script is now geared towards SRR list.
                    # If sra_ids can be SRX, you'd need to use sra-toolkit (e.g. esearch/efetch)
                    # to convert SRX to SRR list before writing to this TSV, or ensure
                    # your `srr_id_input` to the bash script is always SRRs.

                    outfile.write(f"{gsm_id}\t{sra_ids}\n")
                else:
                    print(f"Warning: Could not parse sample line: {sample_line} in entry {i}", file=sys.stderr)
    else:
        print(f"Error: Dataset 'GSE_sample_info_strings' not found in {hdf5_file_path}", file=sys.stderr)

print(f"GSM and SRA ID list written to {output_tsv_path}")