import csv
import numpy
import random
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import cross_val_score

# Create gene interaction sets form the parsed Go Ontology
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

# For each pair of genes from the Collins dataset, create an ontotype feature vector and record the score
# Get iteritems once for efficiency and to make sure items are returned in order
vectors = []
scores = []
collinsData = numpy.load('data/examples/example-genetic-interactions.npy')
num_terms = len(geneSets.keys())
iteritems = geneSets.iteritems()
for i in range(100):
    for j in range(i + 1, 100):
        g1 = "Gene-{}".format(i+1)
        g2 = "Gene-{}".format(j+1)
        vector = [0]*num_terms
        k = 0
        for _, genes in iteritems:
            if g1 in genes:
                vector[k] += 1
            if g2 in genes:
                vector[k] += 1
            k += 1

        vectors.append(vector)
        scores.append(collinsData[i][j])

# Next I shuffle the vectors and scores randomly before doing cross validation
# to eliminate any bias that may exist within the ordering of the data
tmp = list(zip(vectors, scores))
random.shuffle(tmp)
vectors, scores = zip(*tmp)

# Create a RandomForestRegression on the dataset with cross validation
# Output R2 scores for 300 cross validations
rf = RandomForestRegressor(max_depth=10, n_estimators=300)
print cross_val_score(estimator=rf, X=vectors, y=scores, cv=4, scoring="r2")

rf.fit(vectors, scores)
print rf.score(vectors, scores)
