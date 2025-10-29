import h5py

def inspect_gse_xml_dataset(filepath):
    """
    Inspects the 'GSE_XML_strings' dataset in an HDF5 file.

    Args:
        filepath (str): The path to the HDF5 file.
    """
    dataset_name = 'GSE_XML_strings'
    num_strings_to_preview = 100  # How many strings to show a preview of
    chars_to_preview = 500      # How many characters of each string to show

    try:
        with h5py.File(filepath, 'r') as f:
            print(f"Inspecting dataset '{dataset_name}' in HDF5 file: {filepath}\n")

            if dataset_name in f:
                dataset = f[dataset_name]
                print(f"Dataset: {dataset.name}")
                print(f"  Shape: {dataset.shape}")
                print(f"  Data type: {dataset.dtype}")

                if dataset.attrs:
                    print("  Attributes:")
                    for key, val in dataset.attrs.items():
                        print(f"    - {key}: {val}")
                else:
                    print("  No attributes for this dataset.")

                print(f"\n  Previewing the first {num_strings_to_preview} strings (up to {chars_to_preview} characters each):\n")

                # Check if the dataset is not empty
                # if dataset.size > 0 and dataset.ndim > 0:
                print(dataset.shape)
                for i in range(min(num_strings_to_preview, dataset.shape[0])):
                    try:
                        # Read one string element
                        # h5py object dtype datasets often store bytes, so decode if necessary
                        string_data = dataset[i]
                        if isinstance(string_data, bytes):
                            string_data = string_data.decode('utf-8', errors='replace')

                        print(f"--- String {i+1} ---")
                        print(string_data[:chars_to_preview])
                        # if len(string_data) > chars_to_preview:
                        #     print("... [string truncated]\n")
                        # else:
                        print("\n")
                    except Exception as e:
                        print(f"  Error reading or decoding string element {i}: {e}\n")
        #     elif dataset.size == 0:
            #         print("  Dataset is empty.")
            #     else: # scalar dataset or other unexpected shape
            #         try:
            #             string_data = dataset[()]
            #             if isinstance(string_data, bytes):
            #                 string_data = string_data.decode('utf-8', errors='replace')
            #             print(f"--- Scalar String Data ---")
            #             print(string_data[:chars_to_preview])
            #             if len(string_data) > chars_to_preview:
            #                 print("... [string truncated]\n")
            #             else:
            #                 print("\n")
            #         except Exception as e:
            #             print(f"  Error reading or decoding scalar dataset: {e}\n")

            # else:
            #     print(f"Dataset '{dataset_name}' not found in the file.")

    except FileNotFoundError:
        print(f"Error: File not found at {filepath}")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == '__main__':
    # Replace with the path to your HDF5 file
    hdf5_file_path = '/scratch/dhan/order_of_splicing/data/human_expression_profiling_by_HT_seq_gds_metadata.hdf5'
    inspect_gse_xml_dataset(hdf5_file_path)