import pickle
import plotly.graph_objs as go
from plotly.offline import download_plotlyjs, init_notebook_mode, plot
from collections import defaultdict

#read in pickled predictions vs scores
scores = []
for i in range(4):
    scores = scores + pickle.load(open("results{}".format(i), "rb"))

# build hash table of measured scores based on predicted score, as in paper
buckets = defaultdict(list)
for score in scores:
    if score[0] < -0.4:
        buckets[-0.4].append(score[1])
    elif score[0] < -0.32:
        buckets[-0.32].append(score[1])
    elif score[0] < -0.24:
        buckets[-0.24].append(score[1])
    elif score[0] < -0.16:
        buckets[-0.16].append(score[1])
    elif score[0] < -0.08:
        buckets[-0.08].append(score[1])
    elif score[0] < 0:
        buckets[0].append(score[1])
    elif score[0] < 0.08:
        buckets[0.08].append(score[1])
    elif score[0] < 0.16:
        buckets[0.16].append(score[1])
    elif score[0] < 0.24:
        buckets[0.24].append(score[1])
    else:
        buckets["inf"].append(score[1])

# now create box plot for each bucket
boxes = []
for key in sorted(buckets.keys()):
    boxes.append(go.Box(y=buckets[key], name="< {}".format(key)))
layout = go.Layout(
    yaxis=dict(title="Measured Interaction Score", range=[-5, 5]),
    xaxis=dict(title="Predicted Interaction Score")
)

plot(go.Figure(data=boxes, layout=layout))
