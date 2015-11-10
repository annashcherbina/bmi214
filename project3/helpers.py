import math 
def dist():
        return math.sqrt(math.pow(a[0]-b[0],2)+math.pow(a[1]-b[1],2)+math.pow(a[2]-b[2],2))
def main():
        global a
        global b 
        a=[1,3,7]
        b=[1,2,4]
        d=dist()
        print str(d) 
if __name__=="__main__":
        main() 
