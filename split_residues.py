import os

def split_residues(pdb_path, st):
    dirName = pdb_path[:-4]+"_unfolded"
    
    if not os.path.exists(dirName):
        os.mkdir(dirName)
        print("Directory " , dirName ,  " Created ")
    else:    
        print("Directory " , dirName ,  " already exists")
        return -1

    f = open(pdb_path,"r")
    lines = f.readlines()

    index_res = 1
    index_at = 0

    for res in st.get_residues():
        res_f = open(os.path.join(dirName,str(index_res)+".pdb"),"x")
        index_res += 1

        for at in res.get_atoms():
            res_f.write(lines[index_at])
            index_at += 1

    return 1