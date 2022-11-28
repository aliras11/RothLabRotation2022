'''the purpose of this script was to clean the chek2variant_carriers.csv that originally came from the CARRIERS populationbased CHEK2 variant carriers.xlsx file
 '''

import pandas as pd
import re

aa_codes = {'C':'Cys', 'D':'Asp', 'S':'Ser','Q':'Gln',  'K':'Lys',
     'I':'Ile', 'P':'Pro', 'T':'Thr',  'F':'Phe',  'N':'Asn', 
      'G':'Gly',  'H':'His',  'L':'Leu',  'R':'Arg', 'W':'Trp', 
      'A':'Ala', 'V':'Val',  'E':'Glu', 'Y':'Tyr', 'M':'Met'}



chek2_seq_list = []
with open("/Users/alirezarasoulzadeh/Desktop/Roth Lab/aa_seq_chek2.csv") as file:
    next(file)
    for row in file:
       container = row.split(sep=",")
       amino_acid = re.sub(r'"','',container[0])
       number = int(re.sub(r'\D','',container[1]))
       chek2_seq_list.append((number,amino_acid))
chek2_seq_dict = dict(chek2_seq_list)
print(chek2_seq_dict)

#matches the amino acid at the calculated position (nucleotide seq divided by 3) to its three letter amino acid code and returns the appropriately formatted string for matching to the imputed scores file 
def amino_acid_matcher(num):
    aa = chek2_seq_dict[num]
    aa2 = aa_codes[aa]
    return f"p.{aa2}{num}="

def stop_cleaner(string_code):
    p = re.compile(r"p.\w*")
    results = p.search(string_code)
    code = results[0]
    code = re.sub(r"X",'*',code)
    return code

def missense_cleaner(string_code):
    p = re.compile(r"p.\w*")
    results = p.search(string_code)
    code = results[0]
    return code

       
df = pd.read_csv("/Users/alirezarasoulzadeh/Desktop/Roth Lab/chek2variant_carriers.csv",index_col=0)

df['CHEK2 variant_cleaned'] = df['CHEK2 variant'].apply(lambda x: re.sub(r'\D*','',str(x)))

stop_df = df.loc[df["variant type"]=="stop_gained"].copy()
stop_df["hgvs_pro"] =stop_df["CHEK2 variant"].apply(stop_cleaner)
stop_df.to_csv("../nonsense_chek2.csv")


synonymous_df = df.loc[df["variant type"]=="synonymous_variant"].copy()
synonymous_df['aa_position'] =synonymous_df["CHEK2 variant_cleaned"].apply(lambda x: int(x)//3)
synonymous_df["hgvs_pro"] =synonymous_df["aa_position"].apply(amino_acid_matcher)
synonymous_df.to_csv("../synonymous_chek2.csv")


missense_df = df.loc[df["variant type"]=="missense_variant"]
missense_df["hgvs_pro"] =missense_df["CHEK2 variant"].apply(missense_cleaner)
missense_df.to_csv("../missense_chek2.csv")


df.to_csv("chek2_variant_cleaned.csv")