from parameters import *
import sys 


def align_global(seqA,seqB,dx,ex,dy,ey,match_matrix):
    return None



#BUILDS A DICTIONARY OF MATCH SCORES 
def build_match_matrix(pairs):
    match_matrix=dict()
    for line in pairs:
        tokens=line.split(' ')
        token1=tokens[2]
        token2=tokens[3]
        score=int(tokens[4])
        match_matrix[tuple([token1,token2])]=score
        #same score for reverse? Not true, may be a non-symmetric score matrix !
        #match_matrix[tuple([token2,token1])]=score
    return match_matrix

    

#READS IN THE INPUT FILE
def parse_args(data):
    seqA=data[0]
    seqB=data[1]
    local=int(data[2])
    penalties=data[3].split(' ')
    dx=int(penalties[0]) #gap open penalty for A 
    ex=int(penalties[1]) #gap extension penalty for A 
    dy=int(penalties[2]) #gap open penalty for B
    ey=int(penalties[3]) #gap extension penalty for B
    n_alphabet_A=int(data[4])
    alphabet_A=data[5] 
    n_alphabet_B=int(data[6]) 
    alphabet_B=data[7]
    match_matrix=build_match_matrix(data[8::])
    return [seqA,seqB,local,penalties,dx,ex,dy,ey,n_alphabet_A,n_alphabet_B,alphabet_A,alphabet_B,match_matrix]

def main():
    #READ IN THE INPUT FILE 
    data=open(sys.argv[1],'r').read().split('\n')
    while '' in data:
        data.remove('')
    #EXTRACT INPUTS 
    [seqA,seqB,local,penalties,dx,ex,dy,ey,n_alphabet_A,n_alphabet_B,alphabet_A,alphabet_B,match_matrix]=parse_args(data)

    
    #WRITE THE OUTPUT FILE
    outf=open(sys.argv[2],'w')
    

if __name__=="__main__":
    main() 
