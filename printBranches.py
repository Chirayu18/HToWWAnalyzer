import uproot
import argparse

def list_branches(file_path):
    with uproot.open(file_path) as file:
        tree = file["Events"]
        for branch in tree.keys():
            print(branch)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="List all branches from a NanoAOD file.")
    parser.add_argument("file", type=str, help="Path to the NanoAOD ROOT file.")
    args = parser.parse_args()

    list_branches(args.file)
