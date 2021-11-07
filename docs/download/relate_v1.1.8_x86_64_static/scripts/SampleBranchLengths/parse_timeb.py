
import sys
import gzip
import numpy as np

def parse_clues(filename):
    with gzip.open(filename, 'rb') as fp:
        try:
            #parse file
            data = fp.read()
        except OSError:
            with open(filename, 'rb') as fp:
                try:
                    #parse file
                    data = fp.read()
                except OSError:
                    print('Error: Unable to open ' + filename)
                    exit(1)
           
        #get #mutations and #sampled trees per mutation
        filepos = 0
        num_muts, num_sampled_trees_per_mut = np.frombuffer(data[slice(filepos, filepos+8, 1)], dtype = np.int32)
        print(num_muts, num_sampled_trees_per_mut)

        filepos += 8
        #iterate over mutations
        for m in range(0,num_muts):
            bp = np.frombuffer(data[slice(filepos, filepos+4, 1)], dtype = np.int32)
            filepos += 4
            anc, der = np.frombuffer(data[slice(filepos, filepos+2, 1)], dtype = 'c')
            filepos += 2
            daf, n = np.frombuffer(data[slice(filepos, filepos+8, 1)], dtype = np.int32)
            filepos += 8
            print("BP: %d, anc: %s, der %s, DAF: %d, n: %d" % (bp, str(anc), str(der), daf, n))
            
            num_anctimes = 4*(n-daf-1)*num_sampled_trees_per_mut
            anctimes     = np.reshape(np.frombuffer(data[slice(filepos, filepos+num_anctimes, 1)], dtype = np.float32), (num_sampled_trees_per_mut, n-daf-1))
            filepos     += num_anctimes
            print(anctimes)
            
            num_dertimes = 4*(daf-1)*num_sampled_trees_per_mut
            dertimes     = np.reshape(np.frombuffer(data[slice(filepos, filepos+num_dertimes, 1)], dtype = np.float32), (num_sampled_trees_per_mut, daf-1))
            filepos     += num_dertimes
            print(dertimes)

if __name__== "__main__":
    filename = sys.argv[1]
    parse_clues(filename)
