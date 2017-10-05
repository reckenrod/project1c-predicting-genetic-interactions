import csv
from sklearn.ensemble import RandomForestRegressor
import numpy

##Create gene interaction sets form the parsed Go Ontology
geneSets = {}
goMapping = open('data/examples/example-hierarchy-sets.tsv')
goMapping = csv.reader(goMapping, delimiter='\t')
i = 0
for row in goMapping:
    if i not in geneSets:
        geneSets[i] = set()
    for gene in row:
        geneSets[i].add(gene)
    i += 1

##For each pair of genes from the Collins dataset, create an ontotype feature vector and record the score
vectors = []
scores = []
collinsData = numpy.load('data/examples/example-genetic-interactions.npy')
num_terms = len(geneSets.keys())
for i in range(100):
    for j in range(i + 1, 100):
        g1 = "Gene-{}".format(i+1)
        g2 = "Gene-{}".format(j+1)
        vector = [0]*num_terms
        k = 0
        for _, genes in geneSets.iteritems():
            if g1 in genes:
                vector[k] += 1
            if g2 in genes:
                vector[k] += 1
            k += 1

        vectors.append(vector)
        scores.append(collinsData[i][j])

##Train random forest classifier
rf = RandomForestRegressor(max_depth=10, n_estimators=10)
rf.fit(vectors, scores)
print numpy.corrcoef(rf.predict(vectors), scores)
