import matplotlib.pyplot as plt

s = open('DATA/Exactdata.dat', 'r')
T = []
obs = []
data = ((s.read()).split('\n'))

for i in data:
    T.append(float(i.split('\t')[0]))
    obs.append(float(i.split('\t')[1]))

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(T, obs)
