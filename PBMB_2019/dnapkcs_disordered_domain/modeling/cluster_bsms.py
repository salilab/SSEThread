import numpy
from sklearn.cluster import KMeans
from matplotlib import pyplot as plt

'''
Script that clusters an SSEThread output file

Clusters computed using KMeans.

Script is hard-coded to generate only up to three clusters.  Can be modified to produce more.

'''

input_file = "bsms_10000.dat"
plot_cluster_matrix=False
n_clusters = 2
cluster_png_file='pkcs_clustermatrix.png'


def calculate_threading_matrix(seq_length, models, se_lengths):
    # Get % of models with that residue in SE sequence position
    n_se_res = numpy.sum(se_lengths)
    n_se = len(se_lengths)
    thread_matrix = numpy.zeros((n_se_res, seq_length))
    for i in models:
        offset=0
        for s in range(n_se):
            resis = range(int(i[s]), int(i[s]) + se_lengths[s],1)
            for r in range(se_lengths[s]):
                if resis[r] != 0:
                    thread_matrix[offset+r][resis[r]-1]+=1.0

            offset+=se_lengths[s]

    return thread_matrix# / (1.0*len(models))


data = numpy.genfromtxt(input_file, delimiter=" ")

w = data[0:10000,0:3]
km = KMeans(n_clusters).fit(w)

clust0 = numpy.where(km.labels_==0)[0]
clust1 = numpy.where(km.labels_==1)[0]
clust2 = numpy.where(km.labels_==2)[0]

if len(clust0 > 0):
    f = open("cluster0_ids.dat", "w")
    for i in clust0:
        f.write(str(i)+" \n")

if len(clust1 > 0):
    f = open("cluster1_ids.dat", "w")
    for i in clust1:
        f.write(str(i)+" \n")

if len(clust2 > 0):
    f = open("cluster2_ids.dat", "w")
    for i in clust2:
        f.write(str(i)+" \n")

c0 = []
c1 = []
c2 = []

for i in range(len(w)):
    if km.labels_[i] == 0:
        c0.append(data[i])
    elif km.labels_[i] == 1:
        c1.append(data[i])
    elif km.labels_[i] == 2:
        c2.append(data[i])


'''
print("C0-1:", numpy.mean(numpy.array(c0)[:,0]), numpy.std(numpy.array(c0)[:,0]))
print("C0-2:", numpy.mean(numpy.array(c0)[:,1]), numpy.std(numpy.array(c0)[:,1]))
print("C0-3:", numpy.mean(numpy.array(c0)[:,2]), numpy.std(numpy.array(c0)[:,2]))

print("C1-1:", numpy.mean(numpy.array(c1)[:,0]), numpy.std(numpy.array(c1)[:,0]))
print("C1-2:", numpy.mean(numpy.array(c1)[:,1]), numpy.std(numpy.array(c1)[:,1]))
print("C1-3:", numpy.mean(numpy.array(c1)[:,2]), numpy.std(numpy.array(c1)[:,2]))
'''



# Exit if we are not plotting the cluster matrix
if not plot_cluster_matrix:
    exit()

se_lengths = [14,25,10]

tm0 = calculate_threading_matrix(4128, c0, se_lengths)
tm1 = calculate_threading_matrix(4128, c1, se_lengths)
tm2 = calculate_threading_matrix(4128, c2, se_lengths)
#tm = calculate_threading_matrix(seq_length, bsms, se_lengths)

#print numpy.sum(tm[0]), tm[0]

fig, ax=plt.subplots(3,1)
ax[0].imshow(tm0, cmap="Greys", interpolation="none")
ax[0].set_yticklabels([])
ax[0].set_xlabel("Residue Number")
ax[0].set_xlim([2550,2800])

ax[1].imshow(tm1, cmap="Greys", interpolation="none")
ax[1].set_yticklabels([])
ax[1].set_xlabel("Residue Number")
ax[1].set_xlim([2550,2800])

ax[2].imshow(tm2, cmap="Greys", interpolation="none")
ax[2].set_yticklabels([])
ax[2].set_xlabel("Residue Number")
ax[2].set_xlim([2550,2800])


plt.savefig(cluster_png_file, dpi=900)

plt.show()
