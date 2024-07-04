import sys
import pandas as pd


PREFIX = "@"

def main():
    fname = sys.argv[1]
    
    cols_to_keep = ['Antisurvival','Prosurvival']
    
    df = pd.read_csv(fname)
    df.set_index('ID', inplace=True)
    
    columns = df.columns.tolist()
    for c in cols_to_keep:
        columns.remove(c)
        
    for idx in df.index:
        genotype = PREFIX
        genotype += "".join(df.loc[idx,columns].values.astype(str))
        prosurvival = df.Prosurvival[idx]
        antisurvival = df.Antisurvival[idx]
        
        print "%s,%.2f,%.2f" % (genotype,prosurvival,antisurvival)


main()
