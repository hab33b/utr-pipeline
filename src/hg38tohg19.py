import pandas as pd
from pyliftover import LiftOver
import glob


def either(c):
    return '[%s%s]' % (c.lower(), c.upper()) if c.isalpha() else c


def pool_samples(study):
    """
    maf files can not have a title header

    If multiple maf files come from the study, they will be combined into 1 dataframe
    In this case, combining 53 samples from the common cSCC cohort I processed with ngs-pipeline
    and combine these samples with another cohort of 19 primary cSCCs that later metastasized.
    :return: combined maf file
    """
    study_reg = ''.join(either(char) for char in study)     # case insensitivity
    maf_files = glob.glob("../data/mafs-hg38/*" + study_reg + "*")  # finds all studies in directory

    print(maf_files)
    dfs = [pd.read_csv(file, sep="\t", low_memory=False) for file in maf_files]
    dfs = pd.concat(dfs, axis=0)
    print("num of duplicates: ", dfs.duplicated().sum(), "\n", dfs[dfs.duplicated()])
    return dfs.drop_duplicates()


def liftover(maf):
    lo = LiftOver("hg38", "hg19")

    # LIFTOVER
    # 0-index the start_position for LiftOver
    maf["Start_Position"] = maf["Start_Position"].astype("Int64")
    maf["End_Position"] = maf["End_Position"].astype("Int64")

    maf["Start_Position"] = maf["Start_Position"].apply(lambda x: x - 1)
    if "chr" not in maf["Chromosome"].iloc[0]:
        maf["Chromosome"] = maf["Chromosome"].apply(lambda x: "chr" + x)

    # convert positions from hg38 to hg19 using liftover
    maf["Start_Position"] = maf[["Chromosome", "Start_Position", "End_Position"]].apply(
        lambda x: lo.convert_coordinate(x[0], x[1])[0][1] if lo.convert_coordinate(x[0], x[1]) else None, axis=1)
    maf["End_Position"] = maf[["Chromosome", "Start_Position", "End_Position"]].apply(
        lambda x: lo.convert_coordinate(x[0], x[2])[0][1] if lo.convert_coordinate(x[0], x[2]) else None, axis=1)

    # delete failed liftover conversions
    print("LiftOver Conversions failed: ",
          (maf["Start_Position"].isna() | maf["End_Position"].isna()).sum())
    maf = maf[(maf["Start_Position"].isna() | maf["End_Position"].isna()) == False]

    # return position to 1-index
    maf["Start_Position"] = maf["Start_Position"].astype("Int64")
    maf["End_Position"] = maf["End_Position"].astype("Int64")
    maf["Start_Position"] = maf["Start_Position"].apply(lambda x: x + 1)
    return maf


def main():
    # study = input("Enter the TCGA code (cesc, cscc, hnsc, lusc, skcm): ")
    # maf = pool_samples(study)
    # maf = liftover(maf)
    # maf.to_csv("../../data/mafs-hg19/" + study.upper() + "-hg19.maf", sep="\t", index=False)

    # one-maf file
    maf = pd.read_csv("../../data/mafs/mafs-hg38*", sep="\t", low_memory=False)
    maf = liftover(maf)
    maf.to_csv("../../data/mafs/mafs-hg19.maf", sep="\t", index=False)


if __name__ == "__main__":
    main()
