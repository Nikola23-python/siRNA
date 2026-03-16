import pandas as pd
import RNA
import re
from tqdm import tqdm
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq

class FullSiRNAFeatureExtractor:
    def gc_content(self, seq):
        return (seq.count("G")+seq.count("C"))/len(seq)*100
    def fold_energy(self, seq):
        structure, mfe = RNA.fold(seq)
        return mfe
    def duplex_energy(self, s, a):
        d = RNA.duplexfold(s, a)
        return d.energy
    def calc_tm_nn(self, sequence):
        seq_obj = Seq(sequence.replace('U', 'T'))
        tm = mt.Tm_NN(seq_obj, nn_table=mt.RNA_NN1)
        return tm
    def extract(self, s, a):
        s = s.upper().replace("T","U")
        a = a.upper().replace("T","U")
        if not s.endswith("UU"):
            s = s + "UU"
        if not a.endswith("UU"):
            a = a + "UU"
        f = {}

        f["length_19_25"] = 1 if 19 <= len(s) <= 25 else 0 # Индекс 18

        gc = self.gc_content(s)
        f["gc_36_52"] = 1 if 36 <= gc <= 52 else 0 # Индекс 11

        f["no_G_pos13_sense"] = 1 if s[12]!="G" else 0 # Индекс 1

        f["no_GC_pos19_sense"] = 1 if s[18] not in "GC" else 0 # Индекс 2

        f["AU_pos10_sense"] = 1 if s[9] in "AU" else 0 # Индекс 3,34

        f["A_pos3_sense"] = 1 if s[2]=="A" else 0 # Индекс 4

        f["A_pos19_sense"] = 1 if s[18]=="A" else 0 # Индекс 4

        f["GC_pos1_sense"] = 1 if s[0] in "GC" else 0 # Индекс 31

        f["A_pos6_antisense"] = 1 if a[5]=="A" else 0 # Индекс 5

        f["AU_pos1_antisense"] = 1 if a[0] in "AU" else 0 # Индекс 35,42

        f["U_pos2_7_11"] = 1 if (a[1]=="U" and a[6]=="U" and a[10]=="U") else 0 # Индекс 36

        f["C_pos19_antisense"] = 1 if a[18]=="C" else 0 # Индекс 37

        region = a[0:3]
        f["AU_1_3_high"] = 1 if (region.count("A")+region.count("U"))>=2 else 0 # Индекс 38

        f["U_pos13_14"] = 1 if (a[12]=="U" and a[13]=="U") else 0 # Индекс 39

        f["no_A_17_19"] = 1 if "A" not in a[16:19] else 0 # Индекс 40

        f["no_G_14_18"] = 1 if (a[13]!="G" and a[17]!="G") else 0 # Индекс 43

        f["no_U_pos1_sense"] = 1 if s[0] != "U" else 0 # Индекс 32

        first_third = a[:len(a)//3]
        au = first_third.count("A")+first_third.count("U")
        f["AU_first_third_5"] = 1 if au>=5 else 0 # Индекс 29

        region = a[4:19]
        f["AU_5_19_gt6"] = 1 if (region.count("A")+region.count("U"))>6 else 0 # Индекс 33

        antisense_2_7 = a[1:7]
        antisense_8_18 = a[7:18]
        gc1 = self.gc_content(antisense_2_7)
        gc2 = self.gc_content(antisense_8_18)
        tolerance_2_7 = 3
        tolerance_8_18 = 3
        gc_profile_optimal = (abs(gc1 - 19) <= tolerance_2_7) and (abs(gc2 - 52) <= tolerance_8_18)
        f["gc_profile"] = 1 if gc_profile_optimal else 0 # Индекс 10

        f["no_GC_stretch"] = 1 if not re.search(r"[GC]{10,}", s) else 0 # Индекс 30

        f["repeat_limit"] = 1 if s.count("GC")<3 and s.count("AT")<4 else 0 # Индекс 9

        f["motif_UCU_UCCG"] = 1 if ("UCU" in a or "UCCG" in a) else 0 # Индекс 44

        f["no_bad_motif"] = 1 if not any(x in a for x in ["ACGA","GCC","GUGG"]) else 0 # Индекс 45

        f["TT_overhang"] = 1 if s.endswith("UU") else 0 # Индекс 8

        f["duplex_dG"] = self.duplex_energy(s, a)
        f["duplex_dG_point"] = 1 if -35.0 <= self.duplex_energy(s,a) <= -27.0 else 0 # Индекс 46

        f["fold_sense"] = self.fold_energy(s)
        f["fold_sense_point"] = 1 if self.fold_energy(s) > -1.5 else 0 # Индексы 19

        f["fold_antisense"] = self.fold_energy(a)
        f["fold_antisense_point"] = 1 if self.fold_energy(a) > -1.5 else 0 # Индексы 19

        seed_seq = a[1:8]
        tm_value = self.calc_tm_nn(seed_seq)
        f["seed_tm_val"] = tm_value
        f["seed_tm__point"] = 1 if tm_value < 21.5 else 0 # Индекс 20, 21

        valley = s[8:14]
        f["energy_valley"] = self.fold_energy(valley)
        f["energy_valley_point"] = 1 if self.fold_energy(valley) > 0 else 0 # Индекс 17

        dg1 = self.fold_energy(s[0:1])
        dg18 = self.fold_energy(s[17:18])
        f["dg_difference"] = abs(dg1-dg18) # Индекс 41

        f["TOTAL_SCORE"] = sum(v for v in f.values() if isinstance(v,int))

        return f

    def process_dataframe(self, df):
        rows=[]
        for row in tqdm(df.itertuples(index=False), total=len(df)):
            entry=row._asdict()
            try:
                feats=self.extract(
                    entry["sense_sequence"],
                    entry["antisense_sequence"]
                )
                entry.update(feats)
            except:
                entry["TOTAL_SCORE"]=0
            rows.append(entry)
        return pd.DataFrame(rows)

if __name__=="__main__":
    df=pd.read_csv("sirna_merged_all.csv")
    extractor=FullSiRNAFeatureExtractor()
    result=extractor.process_dataframe(df)
    result.to_csv("sirna_full_features.csv",index=False)