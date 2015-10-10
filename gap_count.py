data=open('alignment3.output','r').read().split('\n')
counts=[]
for line in data:
    c=float(line.count('_')) 
    counts.append(c)
meanval=sum(counts)/len(counts)
print str(meanval)

