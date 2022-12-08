from CosmicSqliteapi import *
import re

cosmicdata = CosmicDB("/Users/alirezarasoulzadeh/Desktop/CosmicMutantExport.db")
cosmicdata.open_db()

chek2_data = cosmicdata.gene_query("BRCA1")

#get rid of unnecessary transcript annotation after gene name 
for index in range(len(chek2_data.iloc[:,0])):
    chek2_data.iloc[index,0] = "BRCA1"
    

#get rid of unnecessary annotation in HGVSP column
chek2_data['HGVSP'] = chek2_data['HGVSP'].apply(lambda x: re.sub(r'\w*\.\d:','',str(x)))
print("hello")
chek2_data.to_csv("brca1data_cosmic.csv",index=False)

