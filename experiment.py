import csv
import numpy
import random
import pickle
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import cross_val_score
from scipy.stats import pearsonr

# Create gene interaction sets form the parsed Go Ontology
geneSets = {}
goMapping = open('data/mmc5-gene-term.csv')
goMapping = csv.reader(goMapping, delimiter=',')
for row in goMapping:
    if row[1] not in geneSets:
        geneSets[row[1]] = set()
    geneSets[row[1]].add(row[0])

# For each pair of genes from the Collins dataset,
# create an ontotype feature vector and record the scores
vectors = []
scores = []
collinsData = open('data/collins-sc-emap-gis.tsv')
next(collinsData)
collinsData = csv.reader(collinsData, delimiter='\t')
num_terms = len(geneSets.keys())
for row in collinsData:
    vector = [0]*num_terms
    i = 0
    for _, genes in geneSets.iteritems():
        if row[0] in genes:
            vector[i] += 1
        if row[1] in genes:
            vector[i] += 1
        i += 1
    vectors.append(vector)
    scores.append(float(row[2]))

# Next I shuffle the vectors and scores randomly before doing cross validation
# to eliminate any bias that may exist within the ordering of the data.
# I also split the data into 4 pieces for cross validation during this step.
zipped = list(zip(vectors, scores))
random.shuffle(zipped)
vectors = []
scores = []
z1 = zipped[:len(zipped)/2]
z2 = zipped[len(zipped)/2:]
for z in [z1, z2]:
    v, s = zip(*z[:len(z)/2])
    vectors.append(list(v))
    scores.append(list(s))
    v, s = zip(*z[len(z)/2:])
    vectors.append(list(v))
    scores.append(list(s))

# Create a RandomForestRegression on the dataset with cross validation
# Output pearson scores for 4 cross validations on 10 trees
# Write predictions and scores to file for box plot creation
rf = RandomForestRegressor(n_estimators=10)
for i in range(4):
    trainV = list(vectors)
    trainS = list(scores)
    trainV.pop(i)
    trainS.pop(i)
    rf.fit([item for sublist in trainV for item in sublist], [item for sublist in trainS for item in sublist])
    predictions = rf.predict(vectors[i])
    z = zip(predictions, scores[i])
    pickle.dump(z, open("results{}".format(i), "wb"))
    print(pearsonr(predictions, scores[i]))
