'''
Created on Nov 25, 2013
@author: Ambarish Hazarnis
Run as> python HW5.py
'''

from Sequence import *

# Below function finds- 
#    1. Maximum score for global alignement
#    2. Global Sequence Alignment
#    3. Similarity Measure
   
def allignSequences(S,T):
    S='-'+S
    T='-'+T
    n=len(S)
    m=len(T)
    Seq,V=[],[]  
   
    #Defining the matrices
    for i in range(n):
        V.append([0.0]*m)
        lists=[]
        for j in range(m):
            lists.append(Sequence())
        Seq.append(lists)
    
  
    #Initializing the matrices
    for i in range(n):
        V[i][0]=i*-1
        Seq[i][0].s=S[1:i+1]
        Seq[i][0].t='-'*i 
        
    for j in range(m):
        V[0][j]=j*-1
        Seq[0][j].s='-'*j       
        Seq[0][j].t=T[1:j+1]         
    
    #Dynamic programming approach    
    for i in range(1,n):
        for j in range(1,m):
            ch = getMax(i,j,S,T,Seq,V)
            putMax(ch,i,j,S,T,Seq,V)
            
    #Return result
    return (V[n-1][m-1],Seq[n-1][m-1],getMeasure(Seq[n-1][m-1]))

# Below function finds the position of maximum gain                                                           
def getMax(i,j,S,T,Seq,V):
 
    if S[i] == T[j]:
        a=V[i-1][j-1] + 2
    else:
        a=V[i-1][j-1] - 1
    
    b=V[i-1][j] - 1
    c=V[i][j-1] - 1    
         
    if a >= max(b,c):
        return 1
    elif b >= c:
        return 2
    else:
        return 3

#Below function updates the scores and sequences    
def putMax(ch,i,j,S,T,Seq,V):
    if ch==1:
        if S[i] == T[j]:
            V[i][j]=V[i-1][j-1] + 2
        else:
            V[i][j]=V[i-1][j-1] - 1
        Seq[i][j].s = Seq[i-1][j-1].s + S[i]
        Seq[i][j].t = Seq[i-1][j-1].t + T[j]
    elif ch==2:
        V[i][j]=V[i-1][j] - 1
        Seq[i][j].s = Seq[i-1][j].s + S[i]
        Seq[i][j].t = Seq[i-1][j].t + '-'
    else:
        V[i][j]=V[i][j-1] - 1
        Seq[i][j].s = Seq[i][j-1].s + '-'
        Seq[i][j].t = Seq[i][j-1].t + T[j]
    
    
# Below function finds the similarity measure in the sequence 
def getMeasure(sequence):
    length=len(sequence.s)
    count=0.0
    for i in range(length):
        if sequence.s[i] == sequence.t[i]:
            count+=1
    return count/length


def main():    
    (score,sequence,similarity)=allignSequences("CTGAT", "GTCGA")   
    print "Similarity: " + str(similarity)
    print "Optimal Alignment Score: " + str(score)
    print sequence  
    

main()    