import numpy as np
import csv
import sys
import collections
import itertools as it

#print sys.argv[1]
#print sys.argv[2]              
path=sys.argv[1]

#hardcoding At&T dataset dimensions for now
#dim=180
#dim=5

#e=2
#r=2
e=int(sys.argv[2])
r=float(sys.argv[3])
print "expansion and inflation parameters are"
print e,r

#reading dataset
#reader = csv.reader(open('/Users/saisrinath/Projects/Canopy_Projects/MarkovClustering/new_att.txt'),delimiter=" ")
#reader = csv.reader(open('/Users/saisrinath/Projects/Canopy_Projects/MarkovClustering/new_collaboration.txt'),delimiter=" ")
#reader = csv.reader(open('/Users/saisrinath/Projects/Canopy_Projects/MarkovClustering/new_yeast.txt'),delimiter=" ")

reader = csv.reader(open(path),delimiter=" ")

l=list(reader)

vertices=set()

#creating association matrix

for k,v in l:
    
    vertices.add(k)
    vertices.add(v)

#storing the number of vertices
dim=len(vertices)
print "dim is"
print dim

init = np.zeros((dim,dim),dtype=np.float64) 

indexmap={}     

#using dictionary to map vertex to index
value=0
for i in vertices:
    indexmap[i]=value
    
    value=value+1
    
    
 
#filling the adjancency matrix
for k,v in l:
       
    init[indexmap[k],indexmap[v]]=1
    init[indexmap[v],indexmap[k]]=1              
#using modarr as init without self loops for modularity 
modarr=init
#adding self loops
np.fill_diagonal(init,1) 
#finding sum of columns   
col_sum=init.sum(axis=0,dtype=np.float64)



#normalizing
transition=np.zeros((dim,dim),dtype=np.float64)
#col_sum[col_sum==0]=1
transition=np.divide(init,col_sum)

#transition=normalize(init,col_sum,len)        
#print "transition matrix is"
#print transition   
#expansion


previous=transition


ct=0


while True:
    transition=np.linalg.matrix_power(transition,e) 
    #print transition
    
    #inflation
    
    transition=np.power(transition,r)
    
    
    #normalizing
    
    
    trans_sum=transition.sum(axis=0,dtype=np.float64)    
    
    transition=np.divide(transition,trans_sum)
           
         
    #pruning
    transition=transition.round(4)    
    
    
    ct=ct+1
    print "Iteration" +str(ct)
    
    if (previous==transition).all():
        print "It has converged" 
        break
    #print previous
    previous=transition
    
                                                    
#end of iteration                                                                                                                                              
print "markov matrix is"
print transition


#adding values from the final matrix to a list of clusters
clusters=[]
for t in transition:
    #print i
    indices = [i for i, x in enumerate(t) if x == 1]
    #print indices
    if indices:
        #print indices
        
        minilist=list()
        for valu in indices:
            key = (key for key,value in indexmap.items() if value==valu).next()
            #print key
            minilist.append(key)
            #print minilist
        clusters.append(minilist)
noc=0;            
for c in clusters:
    
      
    noc=noc+1
    print c
   
print "number of clusters are:"
print noc
            

#writing to file       
#f=open('/Users/saisrinath/Projects/Canopy_Projects/MarkovClustering/new_yeast.clu','w')
#removing .txt extension from the path string
path=path[:-4]
f=open(path+'.clu','w')
f.write("*Partition PartitionName\n")
f.write("*Vertices ")
f.write(str(len(vertices))+'\n')

sor=list()
for i in indexmap:
    sor.append(int(i))
sor=sorted(sor)
  
fd={}
for i in indexmap:
    label=1
    for c in clusters:
        
        if i in c:
            fd[int(i)]=label
            #f.write(str(label)+'\n')
           
        label=label+1

for k in sorted(fd):
    #print k,fd[k]
    f.write(str(fd[k])+'\n')
f.close()

#modularity


intra=0
inter=0
np.fill_diagonal(init,0) 
for i in range(dim):
    for j in range(dim):
            
            if init[i,j]==1:               
                key1 = (key for key,value in indexmap.items() if value==i).next()
                
                key2 = (key for key,value in indexmap.items() if value==j).next()
                #print key1,key2
                
                for c in clusters:
                    if key1 in c:
                        if key2 in c:
                            intra=intra+1
                        else:
                            inter=inter+1
                                    
                
print "modularity is:"
print (intra-inter)/2 
                     
                      
        
        
   



  