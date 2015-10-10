import sys 
import itertools 
sofar=[]


def trace_back_local(alignedA,alignedB,pointer_dict,nextkey,seqA,seqB,outf):
    #append to aligned sequences!
    m_index=nextkey[0]
    x_index=nextkey[1]
    y_index=nextkey[2]
    scoreval=nextkey[3]
    if (scoreval==0) and (m_index=="M"):
        dontadd=True
    else: 
        if m_index=="M":
            alignedA=alignedA+seqA[x_index]
            alignedB=alignedB+seqB[y_index]
        elif m_index=="Ix":
            alignedA=alignedA+seqA[x_index]
            alignedB=alignedB+'_'
        elif m_index=="Iy":
            alignedA=alignedA+'_'
            alignedB=alignedB+seqB[y_index]
    if (nextkey not in pointer_dict) or ((scoreval==0) and (m_index=="M")):
        #we're done, break out of the recursion!
        #Trim off the end gaps from sequences and score list
       
        strip_front=0
        strip_rear=0 
        for i in range(len(alignedA)):
            if alignedA[i]=="_":
                strip_front+=1
            elif alignedB[i]=="_":
                strip_front+=1
            else:
                #no more gaps, stop!
                break
        #reverse and check for leading gaps from the other direction        
        alignedA=alignedA[::-1]
        alignedB=alignedB[::-1]
        for i in range(len(alignedA)):
            if alignedA[i]=="_":
                strip_rear+=1
            elif alignedB[i]=="_":
                strip_rear+=1
            else:
                #no more gaps, stop!
                break
        start_index=strip_rear
        end_index=len(alignedA)-strip_front
        alignedA=alignedA[start_index:end_index]
        alignedB=alignedB[start_index:end_index]
        result=tuple([alignedA,alignedB])
        if result not in sofar: 
            outf.write('\n\n'+alignedA+'\n'+alignedB)
            sofar.append(result) 
        return None 
    else:
        #add to the alignd sequences
        steps=pointer_dict[nextkey]
        for s in steps: 
            trace_back_local(alignedA,alignedB,pointer_dict,s,seqA,seqB,outf)


def trace_back(alignedA,alignedB,pointer_dict,nextkey,seqA,seqB,outf):
    #append to aligned sequences!
    m_index=nextkey[0]
    x_index=nextkey[1]
    y_index=nextkey[2]
    
    if m_index=="M":
        alignedA=alignedA+seqA[x_index]
        alignedB=alignedB+seqB[y_index]
    elif m_index=="Ix":
        alignedA=alignedA+seqA[x_index]
        alignedB=alignedB+'_'
    elif m_index=="Iy":
        alignedA=alignedA+'_'
        alignedB=alignedB+seqB[y_index]
    nextkey=nextkey[0:3] 
    if nextkey not in pointer_dict:
        #we're done, break out of the recursion!
        #Trim off the end gaps from sequences and score list
       
        strip_front=0
        strip_rear=0 
        for i in range(len(alignedA)):
            if alignedA[i]=="_":
                strip_front+=1
            elif alignedB[i]=="_":
                strip_front+=1
            else:
                #no more gaps, stop!
                break
        #reverse and check for leading gaps from the other direction        
        alignedA=alignedA[::-1]
        alignedB=alignedB[::-1]
        for i in range(len(alignedA)):
            if alignedA[i]=="_":
                strip_rear+=1
            elif alignedB[i]=="_":
                strip_rear+=1
            else:
                #no more gaps, stop!
                break
        start_index=strip_rear
        end_index=len(alignedA)-strip_front
        alignedA=alignedA[start_index:end_index]
        alignedB=alignedB[start_index:end_index]
        result=tuple([alignedA,alignedB])
        if result not in sofar: 
            outf.write('\n\n'+alignedA+'\n'+alignedB)
            sofar.append(result) 
        return None 
    else:
        #add to the alignd sequences
        steps=pointer_dict[nextkey]
        for s in steps: 
            trace_back(alignedA,alignedB,pointer_dict,s,seqA,seqB,outf)
    
#Local Smith Waterman alignment 
def align_local(seqA,seqB,dx,ex,dy,ey,match_matrix,outf):
    #CREATE THE SCORE MATRIX
    nA=len(seqA)
    nB=len(seqB) 

    M=dict()
    Ix=dict()
    Iy=dict()
    
    p=dict() #pointer matrix 
    
    
    for x in range(nA+1): #add one to account for gap position 
        M[x]=dict()
        p[x]=dict()
        Ix[x]=dict()
        Iy[x]=dict() 
        for y in range(nB+1):
            M[x][y]=None
            p[x][y]=None
            Ix[x][y]=None
            Iy[x][y]=None
            
    #TAKE CARE OF BOUNDARY CONDITIONS
    for x in range(nA+1):
        for y in range(nB+1):
            M[x][y]=0
            Ix[x][y]=0
            Iy[x][y]=0
    pointers=dict()
    
    for x in range(1,nA+1):
        for y in range(1,nB+1): 
            Sx=match_matrix[tuple([seqA[x-1],seqB[y-1]])]
            M_match=round(M[x-1][y-1]+Sx,10)
            M_gapy=round(Ix[x-1][y-1]+Sx,10)
            M_gapx=round(Iy[x-1][y-1]+Sx,10)
            options=[M_match,M_gapy,M_gapx,0]
            score=max(options)
            M[x][y]=score
            slot=tuple(['M',x,y,score])
            
            if M[x][y]==M_match:
                keyval=tuple(['M',x-1,y-1,M[x-1][y-1]])
                if slot in pointers:
                    pointers[slot].append(keyval) 
                else:
                    pointers[slot]=[keyval]
            if M[x][y]==M_gapy:
                keyval=tuple(['Ix',x-1,y-1,Ix[x-1][y-1]])
                if slot in pointers:
                    pointers[slot].append(keyval)
                else:
                    pointers[slot]=[keyval]

            if M[x][y]==M_gapx:
                keyval=tuple(['Iy',x-1,y-1,Iy[x-1][y-1]])
                if slot in pointers:
                    pointers[slot].append(keyval)
                else:
                    pointers[slot]=[keyval]

            Ix_open_gap=round(M[x-1][y]-dy,10)
            Ix_extend_gap=round(Ix[x-1][y]-ey,10)
            options_Ix=[Ix_open_gap,Ix_extend_gap]
            score=max(options_Ix) 
            Ix[x][y]=score 
            slot=tuple(['Ix',x,y,score])
            if Ix[x][y]==Ix_open_gap:
                keyval=tuple(['M',x-1,y,M[x-1][y]])
                if slot in pointers:
                    pointers[slot].append(keyval) 
                else:
                    pointers[slot]=[keyval]
            if Ix[x][y]==Ix_extend_gap:
                keyval=tuple(['Ix',x-1,y,Ix[x-1][y]])
                if slot in pointers:
                    pointers[slot].append(keyval)
                else:
                    pointers[slot]=[keyval]


            Iy_open_gap=round(M[x][y-1]-dx,10)
            Iy_extend_gap=round(Iy[x][y-1]-ex,10)
            Iy_options=[Iy_open_gap,Iy_extend_gap]
            score=max(Iy_options) 
            Iy[x][y]=score 
            slot=tuple(['Iy',x,y,score])
            if Iy[x][y]==Iy_open_gap:
                keyval=tuple(['M',x,y-1,M[x][y-1]])
                if slot in pointers:
                    pointers[slot].append(keyval)
                else:
                    pointers[slot]=[keyval]
            if Iy[x][y]==Iy_extend_gap:
                keyval=tuple(['Iy',x,y-1,Iy[x][y-1]])
                if slot in pointers:
                    pointers[slot].append(keyval)
                else:
                    pointers[slot]=[keyval]



    #DO THE TRACEBACK THROUGH THE POINTER MATRIX!
    max_xy=[] 
    max_score=-1  #all scores in M should be > 0, so we initialized max_score at -1 
    for x in range(nA+1):
        for y in range(nB+1):
            if M[x][y] > max_score:
                max_score=M[x][y]
                max_xy=[tuple([x,y])]
            elif M[x][y]==max_score:
                if tuple([x,y]) not in max_xy:
                    max_xy.append(tuple([x,y]))
    if max_score <=0:
        return None #no alignment can be generated!
    outf.write(str(round(max_score,1)))
    for i in range(len(max_xy)): 
        trace_back_local("","",pointers,tuple(['M',max_xy[i][0],max_xy[i][1],max_score]),'_'+seqA,'_'+seqB,outf)
    return None 

#Global Needleman-Wunsch Alignment 
def align_global(seqA,seqB,dx,ex,dy,ey,match_matrix,outf):
    #CREATE THE SCORE MATRIX
    nA=len(seqA)
    nB=len(seqB) 

    M=dict()
    Ix=dict()
    Iy=dict()
    
    p=dict() #pointer matrix 
    
    
    for x in range(nA+1): #add one to account for gap position 
        M[x]=dict()
        p[x]=dict()
        Ix[x]=dict()
        Iy[x]=dict() 
        for y in range(nB+1):
            M[x][y]=None
            p[x][y]=None
            Ix[x][y]=None
            Iy[x][y]=None
            
    #TAKE CARE OF BOUNDARY CONDITIONS
    #Initialize all as -Inf
    for x in range(nA+1):
        for y in range(nB+1):
            M[x][y]=0
            Ix[x][y]=0
            Iy[x][y]=0 
    #RECURSIVELY FILL IN THE REMAINING POSITIONS
    pointers=dict() 
    for x in range(1,nA+1):
        for y in range(1,nB+1):
            
            Sx=match_matrix[tuple([seqA[x-1],seqB[y-1]])]
            M_match=round(M[x-1][y-1]+Sx,10) 
            M_gapy=round(Ix[x-1][y-1]+Sx,10)
            M_gapx=round(Iy[x-1][y-1]+Sx,10) 
            options=[M_match,M_gapy,M_gapx]
            score=max(options) 
            M[x][y]=score 
            slot=tuple(['M',x,y])
            if M[x][y]==M_match:
                keyval=tuple(['M',x-1,y-1,score])
                if slot in pointers:
                    pointers[slot].append(keyval) 
                else:
                    pointers[slot]=[keyval]
            if M[x][y]==M_gapy:
                keyval=tuple(['Ix',x-1,y-1,score])
                if slot in pointers:
                    pointers[slot].append(keyval)
                else:
                    pointers[slot]=[keyval]

            if M[x][y]==M_gapx:
                keyval=tuple(['Iy',x-1,y-1,score])
                if slot in pointers:
                    pointers[slot].append(keyval)
                else:
                    pointers[slot]=[keyval]

            

            Ix_open_gap=round(M[x-1][y]-dy,10)
            Ix_extend_gap=round(Ix[x-1][y]-ey,10)
            options_Ix=[Ix_open_gap,Ix_extend_gap]
            score=max(options_Ix) 
            Ix[x][y]=score 
            slot=tuple(['Ix',x,y])
            if Ix[x][y]==Ix_open_gap:
                keyval=tuple(['M',x-1,y,score])
                if slot in pointers:
                    pointers[slot].append(keyval) 
                else:
                    pointers[slot]=[keyval]
            if Ix[x][y]==Ix_extend_gap:
                keyval=tuple(['Ix',x-1,y,score])
                if slot in pointers:
                    pointers[slot].append(keyval)
                else:
                    pointers[slot]=[keyval]

                    
            Iy_open_gap=round(M[x][y-1]-dx,10)
            Iy_extend_gap=round(Iy[x][y-1]-ex,10)
            Iy_options=[Iy_open_gap,Iy_extend_gap]
            score=max(Iy_options) 
            Iy[x][y]=score 
            slot=tuple(['Iy',x,y])
            if Iy[x][y]==Iy_open_gap:
                keyval=tuple(['M',x,y-1,score])
                if slot in pointers:
                    pointers[slot].append(keyval)
                else:
                    pointers[slot]=[keyval]
            if Iy[x][y]==Iy_extend_gap:
                keyval=tuple(['Iy',x,y-1,score])
                if slot in pointers:
                    pointers[slot].append(keyval)
                else:
                    pointers[slot]=[keyval]
    ###OUTPUT TEXT FOR QUIZ#####
    print "M:"+str(M)
    print "Ix:"+str(Ix)
    print "Iy:"+str(Iy)
    for p in pointers:
        print str(p)+"-->"+str(pointers[p])

    
    #DO THE TRACEBACK THROUGH THE POINTER MATRIX!
    max_xy=[] 
    max_score=float("-inf")  #all scores in M should be > 0, so we initialized max_score at -1
    #Check last row of M 
    for x in range(nA+1):
        if M[x][nB]> max_score:
            max_score=M[x][nB]
            max_xy=[tuple(['M',x,nB,max_score])]
        elif M[x][nB]==max_score:
            if tuple(['M',x,nB,max_score]) not in max_xy:
                max_xy.append(tuple(['M',x,nB,max_score]))

    #Check last column of M                
    for y in range(nB+1):
        if M[nA][y] > max_score:
            max_score=M[nA][y]
            max_xy=[tuple(['M',nA,y,max_score])]
        elif M[nA][y]==max_score:
            if tuple(['M',nA,y,max_score]) not in max_xy:
                max_xy.append(tuple(['M',nA,y,max_score]))

    startpos=max_xy
    outf.write(str(round(max_score,1)))
    
    alignments=[]
    for sp in startpos:
        #empty alignment and score list 
        trace_back("","",pointers,sp,'_'+seqA,'_'+seqB,outf)
    return None 




#BUILDS A DICTIONARY OF MATCH SCORES 
def build_match_matrix(pairs):
    match_matrix=dict()
    for line in pairs:
        tokens=line.split(' ')
        while '' in tokens:
            tokens.remove('') 
        token1=tokens[2]
        token2=tokens[3]
        score=float(tokens[4])
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
    dx=float(penalties[0]) #gap open penalty for A 
    ex=float(penalties[1]) #gap extension penalty for A 
    dy=float(penalties[2]) #gap open penalty for B
    ey=float(penalties[3]) #gap extension penalty for B
    n_alphabet_A=int(data[4])
    alphabet_A=data[5] 
    n_alphabet_B=int(data[6]) 
    alphabet_B=data[7]
    match_matrix=build_match_matrix(data[8::])
    return [seqA,seqB,local,penalties,dx,ex,dy,ey,n_alphabet_A,n_alphabet_B,alphabet_A,alphabet_B,match_matrix]

def main():
    #READ IN THE INPUT FILE 
    data=open(sys.argv[1],'r').read().split('\n')
    outf=open(sys.argv[2],'w')
    
    while '' in data:
        data.remove('')
    #EXTRACT INPUTS 
    [seqA,seqB,local,penalties,dx,ex,dy,ey,n_alphabet_A,n_alphabet_B,alphabet_A,alphabet_B,match_matrix]=parse_args(data)
    if local==0:
        align_global(seqA,seqB,dx,ex,dy,ey,match_matrix,outf)
    else:
        align_local(seqA,seqB,dx,ex,dy,ey,match_matrix,outf)
    

if __name__=="__main__":
    main() 
