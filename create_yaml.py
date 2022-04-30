import yaml
import sys
import pandas as pd

def create_yaml(filename, output):
    with open(output, "w") as yaml_file:
        
        lines = list(pd.unique(pd.read_table(filename)["patient"]))
        
        data = {}
        data["patients"] = {}
        patients = data["patients"]

        for line in lines:
            line = line.strip()
            patients[line] = {}
            patients[line]["normal"] = ""
            patients[line]["tumor"] = ""
            patients[line]["normal"] = "../bams/" + line + ".normal.bam"
            patients[line]["tumor"] = "../bams/" + line + ".tumor.bam"

        yaml.dump(data, yaml_file, default_flow_style=False)


def main():
    create_yaml("../patients.csv", "config.yaml")


if __name__ == '__main__':
    main()
