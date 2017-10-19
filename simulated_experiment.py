import csv
import numpy
import random
import pickle
from sklearn.ensemble import RandomForestRegressor
from scipy.stats import pearsonr

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
rf = RandomForestRegressor(max_depth=10, n_estimators=10)
for i in range(4):
    trainV = list(vectors)
    trainS = list(scores)
    trainV.pop(i)
    trainS.pop(i)
    rf.fit([item for sublist in trainV for item in sublist], [item for sublist in trainS for item in sublist])
    predictions = rf.predict(vectors[i])
    z = zip(predictions, scores[i])
    pickle.dump(z, open('s_results{}'.format(i), 'wb'))
    print(pearsonr(predictions, scores[i]))
