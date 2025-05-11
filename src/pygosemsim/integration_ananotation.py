import sys

sys.path.append(
    "/net/kihara/scratch/etugoluk/2024_generate_pfp_prediction/final_dataset_constrct_with_viruses/testing_memory/pygosemsim"
)
from pygosemsim import annotation


def main():
    if len(sys.argv) != 2:
        print("Usage: python your_script.py <resource_identifier>")
        sys.exit(1)
        print("test1")

    resource_identifier = sys.argv[1]
    print("test2")
    Whole_annot, y = annotation.from_resource(resource_identifier)
    print("test3")


if __name__ == "__main__":
    main()
