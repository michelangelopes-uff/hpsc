import json
import matplotlib
import matplotlib.pyplot as plt

def readJSON(_filename):
    data = json.load(open(_filename))
    return data

res = readJSON("output.json")
plt.plot(res["resultado"])
plt.show()
