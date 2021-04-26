#First exercise
#Get the zscore of the count of "n = 0, l = 1", "n = 0, l = 2" and "n = 0, l = 3" building blocks on the bacilus network.
library(fibrationSymmetries)
library(dplyr)
library(igraph)

bacilus_network = read.table("./Bacilus.txt", sep = " ", header = F)

buildingBlocks = get.building.blocks(raw_edges = bacilus_network, progressBar = F)
buildingBlocks = filter(buildingBlocks, nl == "n = 0, l = 1")
realCount = nrow(buildingBlocks)

buildingBlocks2 = get.building.blocks(raw_edges = bacilus_network, progressBar = F)
buildingBlocks2 = filter(buildingBlocks2, nl == "n = 0, l = 2")
realCount2 = nrow(buildingBlocks2)

buildingBlocks3 = get.building.blocks(raw_edges = bacilus_network, progressBar = F)
buildingBlocks3 = filter(buildingBlocks3, nl == "n = 0, l = 3")
realCount3 = nrow(buildingBlocks3)

graph <- graph_from_edgelist(as.matrix(bacilus_network))

count = NULL
count2 = NULL
count3 = NULL
sampleSize = 100
for(i in 1:sampleSize) {
  outDegrees <- degree(graph = graph, mode = "out", loops = T, normalized = FALSE)
  outDegrees <- outDegrees[sample(vcount(graph))]
  
  inDegrees <- degree(graph = graph, mode = "in", loops = T, normalized = FALSE)
  inDegrees <- inDegrees[sample(vcount(graph))]
  
  newGraph <- sample_degseq(out.deg = outDegrees, in.deg = inDegrees, method = "simple")
  newEdges <- as.data.frame(as_edgelist(newGraph), stringsAsFactors = F)
  newEdges[] <- apply(newEdges, 2, as.character)
  
  generatedBlocks = get.building.blocks(raw_edges = newEdges, progressBar = F)
  generatedBlocks = filter(generatedBlocks, nl == "n = 0, l = 1")
  count = c(count, nrow(generatedBlocks))
  
  generatedBlocks2 = get.building.blocks(raw_edges = newEdges, progressBar = F)
  generatedBlocks2 = filter(generatedBlocks2, nl == "n = 0, l = 2")
  count2 = c(count2, nrow(generatedBlocks2))
  
  generatedBlocks3 = get.building.blocks(raw_edges = newEdges, progressBar = F)
  generatedBlocks3 = filter(generatedBlocks3, nl == "n = 0, l = 3")
  count3 = c(count3, nrow(generatedBlocks3))
}
print(count)
print(count2)
print(count3)

zscore = (realCount - mean(count)) / sd(count)
print("zscore for n=0, l=1:")
print(zscore)

zscore2 = (realCount2 - mean(count2)) / sd(count2)
print("zscore for n=0, l=2:")
print(zscore2)

zscore3 = (realCount3 - mean(count3)) / sd(count3)
print("zscore for n=0, l=3:")
print(zscore3)

#Second exercise
#Plot the modularity of Newman clusters of a scale-free network vs size of the network for sizes 10-1000 with the step of 10.
clusters_to_dataframe <- function(cluster_data) {
  clusters_df = data.frame(Node = character(), ClusterId = numeric(), stringsAsFactors = F)
  
  for(i in 1:length(cluster_data)) {
    newData = data.frame(Node = cluster_data[[i]], ClusterId = i, stringsAsFactors = F)
    clusters_df = rbind(clusters_df, newData)
  }
  
  clusters_df = arrange(clusters_df, Node)
  
  return(clusters_df)
}

sizes = c()
newmanModularities = c()
for(i in 1:100) {
  sizes = append(sizes, i*10)
  graph = sample_pa(i*10, directed = F)
  newman_clusters = clusters_to_dataframe(cluster_edge_betweenness(graph))
  newmanModularity = modularity(graph, newman_clusters$ClusterId)
  newmanModularities = append(newmanModularities, newmanModularity)
}

plot(sizes, newmanModularities, main="Newman Modularity vs Network Size")