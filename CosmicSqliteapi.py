
import sqlite3 as sq
import os.path 
import pandas as pd

class DataBase:
    def __init__(self, path):
        self.path = path
        self.db = None
        self.cursor = None

    def open_db(self):
        if os.path.exists(os.path.abspath(self.path)):
            try:
                self.path = os.path.abspath(self.path)
                self.db = sq.connect(self.path)
                self.cursor = self.db.cursor()
                self.cursor.execute('PRAGMA Foreign_Keys = True')
            except:
                print("Error Connecting to Database...")
                self.cursor = None
                raise
        else:
            print("Database is not present on system")


    def close_db(self):
        try:
            print("Closing database...")
            self.db.commit()
            self.db.close()
        except:
            print("Error closing database")
            raise


    @property
    def tables(self):
        query = '''PRAGMA table_list;'''
        tablesinfolist =  self.cursor.execute(query).fetchall()
        df = pd.DataFrame(data=None, columns=["Table_Name","Number_of_Columns"])
        for tableinfo in tablesinfolist:
            temp_dict = {"Table_Name":tableinfo[1],"Number_of_Columns":tableinfo[3]}
            tempseries = pd.Series(temp_dict)
            df.loc[len(df)] = tempseries
        return df  #returns a pandas dataframe containing the information outlined in temp_dict


    def table_info(self,table):
        query = f'''PRAGMA table_info({table});'''
        colinfolist =  self.cursor.execute(query).fetchall()
        df = pd.DataFrame(data=None, columns=["Column_Index","Column_Name","Column_Datatype"])
        for col_info in colinfolist:
            temp_dict = {"Column_Index":col_info[0],"Column_Name":col_info[1],"Column_Datatype":col_info[2]}
            tempseries = pd.Series(temp_dict)
            df.loc[len(df)] = tempseries
        return df #returns a pandas dataframe containing the information outlined in temp_dict

    def table_head(self,table,count=50):
        query = f"""select * from {table} limit {count};"""
        tableheadlist =  self.cursor.execute(query).fetchall() #fetchall returns a list of tuples
        df = pd.DataFrame(tableheadlist,columns=list(self.table_info(table).iloc[:,1])) #assumes that column names for the table are in the second column of the pandas df returned by calling table_info
        return df




class CosmicDB(DataBase):

    def gene_query(self,gene,table="MutantExport",size_output=100000):
        query = f'''select * from {table} where gene_name like \'{gene}%\''''
        genequerylist = self.cursor.execute(query).fetchmany(size=size_output)
        df = pd.DataFrame(genequerylist,columns=list(self.table_info(table).iloc[:,1]))
        return df

    #there may be a chance that even after query the resulting array is too big so there is a default max size of 100,000 applied, it is possible
    def organ_query(self,tissue,table="MutantExport",size_output=100000):
        query = f'''select * from {table} where Primary_site like \'{tissue}%\''''
        genequerylist = self.cursor.execute(query).fetchmany(size=size_output)
        df = pd.DataFrame(genequerylist,columns=list(self.table_info(table).iloc[:,1]))
        return df

    def organ_gene_query(self,tissue,gene,table="MutantExport",size_output=100000):
        query = f'''select * from {table} where Primary_site like \'{tissue}%\' and gene_name like \'{gene}%\''''
        genequerylist = self.cursor.execute(query).fetchmany(size=size_output)
        df = pd.DataFrame(genequerylist,columns=list(self.table_info(table).iloc[:,1]))
        return df
    
    
    
    






if __name__ == "__main__":
   testdb = CosmicDB("/Users/alirezarasoulzadeh/Desktop/CosmicMutantExport.db")
   testdb.open_db()
   print(testdb.path)
   print(testdb.tables)
   print(testdb.table_info("MutantExport"))
   print(testdb.gene_query("NCOR2"))
   print(testdb.organ_gene_query("liver","NCOR2"))
   testdb.close_db()