
import pandas as pd
from pyliftover import LiftOver


def main():
    variants = pd.read_csv('variants.maf', sep="\t")
    lo = LiftOver("hg38", "hg19")

    # 0-index bed file
    variants["Start_Position"] = variants["Start_Position"].apply(lambda x: x - 1)

    # convert positions from hg38 to hg19
    variants["Start_Position"] = variants[["Chromosome", "Start_Position", "End_Position"]].apply(
        lambda x: lo.convert_coordinate(x[0], x[1])[0][1] if lo.convert_coordinate(x[0], x[1]) else None, axis=1)
    variants["End_Position"] = variants[["Chromosome", "Start_Position", "End_Position"]].apply(
        lambda x: lo.convert_coordinate(x[0], x[2])[0][1] if lo.convert_coordinate(x[0], x[2]) else None, axis=1)

    # subset & rename
    variants = variants[(variants["Start_Position"].isna() | variants["Start_Position"].isna()) == False]
    mut = variants[["Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode"]]
    mut["Start_Position"] = mut["Start_Position"].astype(int)
    mut.columns = ["CHROMOSOME", "POSITION", "REF", "ALT", "SAMPLE"]

    # utr3 & utr5
    utr3 = mut[variants["Variant_Classification"] == "3'Flank"]
    utr5 = mut[variants["Variant_Classification"] == "5'Flank"]

    utr3.to_csv("variants_utr3", sep="\t", index=False)
    utr5.to_csv("variants_utr5", sep="\t", index=False)


if __name__ == "__main__":
    main()
