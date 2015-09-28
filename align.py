import sys 
import itertools 

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
        outf.write('\n\n'+alignedA+'\n'+alignedB)
        return None #collapse the recursion! 
    
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
        Ix[x][0]=float("-inf")
        Iy[x][0]=float("-inf")

    for y in range(nB+1):
        Ix[0][y]=float("-inf")
        Iy[0][y]=float("-inf")
    
    for x in range(nA+1):
        M[x][0]=0
    for y in range(nB+1):
        M[0][y]=0
        
    pointers=dict()
    
    for x in range(1,nA+1):
        for y in range(1,nB+1): 
            Sx=match_matrix[tuple([seqA[x-1],seqB[y-1]])]
            M_match=M[x-1][y-1]+Sx
            M_gapy=Ix[x-1][y-1]+Sx
            M_gapx=Iy[x-1][y-1]+Sx
            options=[M_match,M_gapy,M_gapx,0]
            score=max(options) 
            M[x][y]=score
            if score >0: 
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
            #if value is 0, we stop the traceback! 
            Ix_open_gap=M[x-1][y]-dy
            Ix_extend_gap=Ix[x-1][y]-ey
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


            Iy_open_gap=M[x][y-1]-dx
            Iy_extend_gap=Iy[x][y-1]-ex
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
        trace_back("","",pointers,tuple(['M',max_xy[i][0],max_xy[i][1],max_score]),'_'+seqA,'_'+seqB,outf)
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
        M[x][0]=float("-inf")
        Ix[x][0]=float("-inf")
        Iy[x][0]=float("-inf")

    for y in range(nB+1):
        M[0][y]=float("-inf")
        Ix[0][y]=float("-inf")
        Iy[0][y]=float("-inf")
    
    for x in range(1,nA+1):
        Ix[x][0]=-dy-(x-1)*ey
    for y in range(1,nB+1):
        Iy[0][y]=-dx-(y-1)*ex
        
    M[0][0]=0
    
    #RECURSIVELY FILL IN THE REMAINING POSITIONS
    pointers=dict() 
    for x in range(1,nA+1):
        for y in range(1,nB+1):
            
            Sx=match_matrix[tuple([seqA[x-1],seqB[y-1]])]
            M_match=M[x-1][y-1]+Sx
            M_gapy=Ix[x-1][y-1]+Sx
            M_gapx=Iy[x-1][y-1]+Sx
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

            

            Ix_open_gap=M[x-1][y]-dy
            Ix_extend_gap=Ix[x-1][y]-ey
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

                    
            Iy_open_gap=M[x][y-1]-dx
            Iy_extend_gap=Iy[x][y-1]-ex
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
            
            
    #DO THE TRACEBACK THROUGH THE POINTER MATRIX!
    max_endscore=max([M[nA][nB],Ix[nA][nB],Ix[nA][nB]])
    if max_endscore==M[nA][nB]:
        startpos.append(tuple(['M',nA,nB,max_endscore]))
    if max_endscore==Ix[nA][nB]:
        startpos.append(tuple(['Ix',nA,nB,max_endscore]))
    if max_endscore==Iy[nA][nB]:
        startpos.append(tuple(['Iy',nA,nB,max_endscore]))
    outf.write(str(round(max_endscore,1)))
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
    outf=open(sys.argv[2],'w')
    
    while '' in data:
        data.remove('')
    #EXTRACT INPUTS 
    [seqA,seqB,local,penalties,dx,ex,dy,ey,n_alphabet_A,n_alphabet_B,alphabet_A,alphabet_B,match_matrix]=parse_args(data)

    if local==False:
        align_global(seqA,seqB,dx,ex,dy,ey,match_matrix,outf)
    else:
        align_local(seqA,seqB,dx,ex,dy,ey,match_matrix,outf)
    

if __name__=="__main__":
    main() 
