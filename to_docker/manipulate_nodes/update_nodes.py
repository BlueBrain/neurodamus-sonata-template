import h5py
import sys

def is_bytes(obj):
    if isinstance(obj[0], bytes):
        return True
    else:
        return False

def update_datasets(nodes_file, dataset_value_pairs):
    with h5py.File(nodes_file, "r+") as f:
        # Load the population from the HDF5 file
        populations = list(f["nodes"].keys())
        assert len(populations) == 1
        pop = populations[0]
        print(f"Population: {pop}")

        for dataset_name, new_param_value in dataset_value_pairs:
            # Check if the dataset exists in the HDF5 file
            if dataset_name not in f["nodes"][pop]['0']:
                print(f"Dataset '{dataset_name}' not found in nodes.h5 file")
                continue

            # Get the dataset and its current value
            dataset = f["nodes"][pop]['0'][dataset_name]
            current_value = dataset[:]

            # Determine the data type of the dataset and convert the new value accordingly
            if is_bytes(current_value):
                new_value = new_param_value.encode()
            else:
                new_value = float(new_param_value) if '.' in new_param_value else int(new_param_value)

            # Update the dataset with the new value
            dataset[:] = new_value
            print(f"Updated {dataset_name} from {current_value} to {new_value}")

if __name__ == "__main__":
    nodes_file = sys.argv[1]
    dataset_value_pairs = [(sys.argv[i], sys.argv[i+1]) for i in range(2, len(sys.argv), 2)]

    update_datasets(nodes_file, dataset_value_pairs)