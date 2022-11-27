from CosmicSqliteapi import *
import re

cosmicdata = CosmicDB("/Users/alirezarasoulzadeh/Desktop/CosmicMutantExport.db")
cosmicdata.open_db()

chek2_data = cosmicdata.gene_query("CHEK2")

#get rid of unnecessary transcript annotation after gene name 
for index in range(len(chek2_data.iloc[:,0])):
    chek2_data.iloc[index,0] = "CHEK_2"
    print(chek2_data.iloc[index,0])


chek2_data['HGVSP'] = chek2_data['HGVSP'].apply(lambda x: re.sub(r'\w*\.\d:','',str(x)))

chek2_data.to_csv("../chek2data_cosmic.csv",index=False)


