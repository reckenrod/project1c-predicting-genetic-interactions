import csv
from sklearn.ensemble import RandomForestRegressor
import numpy

##Create gene interaction sets form the parsed Go Ontology
geneSets = {}
goMapping = open('data/mmc5-gene-term.csv')
goMapping = csv.reader(goMapping, delimiter=',')
for row in goMapping:
    if row[1] not in geneSets:
        geneSets[row[1]] = set()
    geneSets[row[1]].add(row[0])

##For each pair of genes from the Collins dataset, create an ontotype feature vector and record the score
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

##Train random forest classifier on whole dataset and find r2 score
rf = RandomForestRegressor(max_depth=4, n_estimators=4)
rf.fit(vectors, scores)
print numpy.corrcoef(rf.predict(vectors), scores)
