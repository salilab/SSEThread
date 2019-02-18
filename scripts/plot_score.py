import json
import glob
from matplotlib import pyplot as plt 
import numpy as np 
import pandas
from pandas.io.json import json_normalize
import math 
psipred = []
score=[]

for outfile in glob.glob("enumerate_multi_1_15000.dat"):
    with open(outfile) as f:
        lst = json.load(f)
        #print type(lst), type(lst['models'])
        #df = pandas.DataFrame(lst['models'])
        all_rest_scores = []
        for m in lst['models']:
            if not  math.isnan(m['score']):
                rest_scores = []
                for k in m['restraints'].keys():
                    rest_scores.append(m['restraints'][k])
                all_rest_scores.append(rest_scores)
                #psipred.append((m['restraints']['SecondaryStructureParsimonyRestraint0']-1.738)*500)
                score.append(m['score'])

                #if m['score'] > 131 and m['score'] < 132:
                #    print m
                #exit()


#print all_rest_scores.shape

df = pandas.DataFrame(all_rest_scores, columns=m['restraints'].keys())

print len(score), "models"

#plt.hist(np.array(score)[~np.isnan(np.array(score))])
#plt.show()

print min(score), max(score)

for f in df.columns:
    print f, min(df[f]), max(df[f])
#print min(psipred), max(psipred)

hist = df.hist()

plt.show()