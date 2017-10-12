import csv
import numpy
import random
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import cross_val_score

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
# to eliminate any bias that may exist within the ordering of the data
tmp = list(zip(vectors, scores))
random.shuffle(tmp)
vectors, scores = zip(*tmp)

# Create a RandomForestRegression on the dataset with cross validation
# Output R2 scores for 10 cross validations
rf = RandomForestRegressor(max_depth=10, n_estimators=10)
print cross_val_score(estimator=rf, X=vectors, y=scores, cv=10, scoring="r2")
