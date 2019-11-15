#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 16:47:36 2019

@author: Ayca
"""
import numpy as np
import math
#%% File Reading 
def fileClean(string):
    s = 'AGCT'
    result = ''.join([i for i in string if i in s])
    return result

def fileRead(filepath):
    p = ""
    with open(filepath) as fp:
       line = fp.readline()
       cnt = 1
       while line:
           line = fp.readline()
           
           cnt += 1
           if line == ">another_sequence":
               break
           p = p+line
    p = p.split(">another_sequence")
    for i in range(2):
        p[i] = fileClean(p[i])
    
    s1 = p[0]
    s2 = p[1]
    return s1,s2

def fileWrite(score,s1, s2):
    count = 60
    begin = 0
    with open('result.txt', 'w') as writer:
        st = "Score = " + str(score)
        writer.write(st)
        line_len = math.ceil(len(s1)/60)
        for line in range(0, line_len):
            writer.write("\n")
            writer.write("my_first_sequence\t")
            writer.write(s1[begin:count])
            writer.write("\nanother_sequence\t")
            writer.write(s2[begin:count])
            begin = count
            if (count+60 > len(s1)):
                count = len(s1)
            else:
                count = count+60
            writer.write("\n")
    writer.close()
#%%
#scoringMatrix = [[2,-3,-3,-3],
#                 [-3,2,-3,-3],
#                 [-3,-3,2,-3],
#                 [-3,-3,-3,2]]

match_score = 2
mismatch_score = -3
#  A    C   G       T 
#A 2   -3   -3      -3 
#C -3   2   -3      -3 
#G -3  -3    2      -3 
#T -3  -3   -3       2

#match=2
#mismatch = -3
def Matrix(s1,s2,gapOpen):
    row = len(s1)+1
    col = len(s2)+1
    matrix = np.zeros((row,col), dtype=int)
    for i in range(1,row):
        matrix[i][0]=matrix[i-1][0] + gapOpen
    for i in range(1,col):
        matrix[0][i]=matrix[0][i-1] + gapOpen
    return matrix

def LocalMatrix(s1,s2):
    row = len(s1)+1
    col = len(s2)+1
    matrix = np.zeros((row,col), dtype=int)
    return matrix
    
def EMatrix(s1,s2,gapOpen,gapExt):
    row = len(s1)+1
    col = len(s2)+1
    e = np.zeros((row,col), dtype=int)
    e[0][0] = -10000
    for i in range(1,row):
        e[i][0] = gapOpen+i*gapExt
    return e

def FMatrix(s1,s2,gapOpen,gapExt):
    row = len(s1)+1
    col = len(s2)+1
    f = np.zeros((row,col), dtype=int)
    f[0][0] = -10000
    for j in range(1,col):
        f[0][j] = gapOpen+j*gapExt
    return f

def GMatrix(s1,s2,gapOpen,gapExt):
    row = len(s1)+1
    col = len(s2)+1
    g = np.zeros((row,col), dtype=int)
    g[0][0] = -10000
    return g

def VMatrix(s1,s2,gapOpen,gapExt):
    row = len(s1)+1
    col = len(s2)+1
    v = np.zeros((row,col), dtype=int)
    v[0][0] = 0
    for i in range(1,row):
        v[i][0] = gapOpen+i*gapExt
    for j in range(1,col):
        v[0][j] = gapOpen+j*gapExt
    return v

def Score(c1,c2):
    if c1 == c2:
        return match_score
    else:
        return mismatch_score
    
#%%
def Align(matrix, s1,s2, gapOpen,i,j):
    align1 = ""
    align2 = ""
    count = 0
    while i >= 1 and j >= 1 and count <=5:
        current_score = matrix[i][j]
        if current_score == matrix[i-1][j-1] + Score(s1[i-1], s2[j-1]):
            align1 = align1+ s1[i-1]
            align2 = align2+ s2[j-1]
            i -=1
            j -=1
        elif current_score == matrix[i][j-1] + gapOpen:
            align1 = align1 + s1[i]
            align2 = align2 + "-"
            j -=1
        elif current_score == matrix[i-1][j] + gapOpen:
            align1 = align1 + "-"
            align2 = align2 + s2[j]
            i -=1
        count += 1
    align1 = align1[::-1]
    align2 = align2[::-1]
    print("align1: {}, align2: {}".format(align1,align2))
    return align1, align2

def AAlign(matrix, s1,s2, gapOpen, gapExt,i,j):
    align1 = ""
    align2 = ""
    mis_match = 0
    count = 0
    while i >= 1 and j >= 1:
        current_score = matrix[i][j]
        if s1[i-1] == s2[j-1]:
            mis_match = 2
        else:
            mis_match = -3
        if current_score == matrix[i-1][j-1] + mis_match:
            align1 = align1+ s1[i-1]
            align2 = align2+ s2[j-1]
            i -=1
            j -=1

        elif (current_score == matrix[i][j-1] + gapExt) or (current_score == matrix[i][j-1] + gapOpen+ gapExt):
            align1 = align1 + s1[i]
            align2 = align2 + "-"
            j -=1
        elif (current_score == matrix[i-1][j] + gapExt) or (current_score == matrix[i-1][j] + gapOpen+ gapExt):
            align1 = align1 + "-"
            align2 = align2 + s2[j]
            i -=1
        count += 1
    align1 = align1[::-1]
    align2 = align2[::-1]
    print("align1: {}, align2: {}".format(align1,align2))
    return align1, align2

def LocalAlign(matrix, s1,s2, gapOpen,i,j):
    align1 = ""
    align2 = ""
    current_score = matrix[i][j]

    while i >= 1 and j >= 1 and current_score >0:
        current_score = matrix[i][j]

        if s1[i-1] == s2[j-1]:
            mis_match = 2
        else:
            mis_match = -3
        print("\ncurrent_score",current_score)
        print("for Diagonal: matrix[i-1][j-1]:{} + mis_match:{} = {}".format(matrix[i-1][j-1],mis_match, matrix[i-1][j-1] + mis_match))
        print("for Up : matrix[i][j-1]:{} + gapOpen:{} = {}".format(matrix[i][j-1],gapOpen, matrix[i-1][j-1] + gapOpen))
        print("for Left : matrix[i-1][j]:{} + gapOpen:{} = {}".format(matrix[i-1][j],gapOpen, matrix[i-1][j-1] + gapOpen))
        if current_score == matrix[i-1][j-1] + mis_match:
            align1 = align1+ s1[i-1]
            align2 = align2+ s2[j-1]
            i -=1
            j -=1
            print("\nDiagonal chosen")
        elif current_score == matrix[i][j-1] + gapOpen:
            align1 = align1 + s1[i]
            align2 = align2 + "-"
            j -=1
            print("\nUp chosen")
        elif current_score == matrix[i-1][j] + gapOpen:
            align1 = align1 + "-"
            align2 = align2 + s2[j]
            i -=1
            print("\nLeftChosen")
    align1 = align1[::-1]
    align2 = align2[::-1]
    print("align1: {}, align2: {}".format(align1,align2))
    return align1, align2

def LocalAAlign(matrix, s1,s2, gapOpen, gapExt,i,j):
    align1 = ""
    align2 = ""
    mis_match = 0
    current_score = matrix[i][j]
    while i >= 1 and j >= 1 and current_score>0:
        current_score = matrix[i][j]
        if s1[i-1] == s2[j-1]:
            mis_match = 2
        else:
            mis_match = -3
        if current_score == matrix[i-1][j-1] + mis_match:
            align1 = align1+ s1[i-1]
            align2 = align2+ s2[j-1]
            i -=1
            j -=1

        elif (current_score == matrix[i][j-1] + gapExt) or (current_score == matrix[i][j-1] + gapOpen+ gapExt):
            align1 = align1 + s1[i]
            align2 = align2 + "-"
            j -=1
        elif (current_score == matrix[i-1][j] + gapExt) or (current_score == matrix[i-1][j] + gapOpen+ gapExt):
            align1 = align1 + "-"
            align2 = align2 + s2[j]
            i -=1
    align1 = align1[::-1]
    align2 = align2[::-1]
    print("align1: {}, align2: {}".format(align1,align2))
    return align1, align2
#%%
def Global(s1, s2,gapOpen):
    row = len(s1)+1
    col = len(s2)+1
    matrix = Matrix(s1, s2,gapOpen)
    for i in range(1,row):
        for j in range(1,col):
            match = 0
            gap_left = 0
            gap_up = 0
            match = matrix[i-1][j-1]+ Score(s1[i-1],s2[j-1])
            gap_left = matrix[i][j-1] +gapOpen
            gap_up = matrix[i-1][j]+gapOpen
            matrix[i][j] = max(match,gap_left,gap_up)
    score = matrix[row-1,col-1]
    print(score)
    
    i = row-1
    j = col-1
    align = Align(matrix, s1,s2,gapOpen,i,j)
    return matrix, align[0], align[1],score

def AGlobal(s1,s2,gapOpen, gapExt):
    row = len(s1)+1
    col = len(s2)+1
    e = EMatrix(s1, s2,gapOpen,gapExt)
    f = FMatrix(s1, s2,gapOpen,gapExt)
    g = GMatrix(s1, s2,gapOpen,gapExt)
    v = VMatrix(s1, s2,gapOpen,gapExt)
    for i in range(1,row):
        for j in range(1,col):
            e[i][j] = max(0,e[i][j-1]+gapExt,v[i][j-1]+gapOpen+gapExt)
            f[i][j] = max(0,f[i-1][j]+gapExt,v[i-1][j]+gapOpen+gapExt)
            g[i][j] = v[i-1][j-1]+ Score(s1[i-1],s2[j-1])
            v[i][j] = max(0,g[i][j],e[i][j],f[i][j])
    score = v[row-1,col-1]
    print(score)
    i = row-1
    j = col-1
    align = AAlign(v, s1,s2, gapOpen, gapExt,i,j)
    return v, align[0], align[1],score


#AGlobal("TCGACCCAAGTAGGGAAAGAATATCAACACAAAGGCTCGAGAAGAGCCACCCCATGAGCCACCGCATCTACCCCGTGCCCCAGCAAATTAAGAATAG","TCGACCCATGTAGGGAAAGCATATCAATTTCACAAAGGCTCGAGAAGAGCCACATGAGCCACCGCATCTACCCCAGCAAATTAAGAAAAG", -5, -2)
def Local(s1,s2, gapOpen):
    row = len(s1)+1
    col = len(s2)+1
    matrix = LocalMatrix(s1, s2)
    score = 0
    largest_position = (0,0)
    for i in range(1,row):
        for j in range(1,col):
            match = 0
            gap_left = 0
            gap_up = 0
            match = matrix[i-1][j-1]+ Score(s1[i-1],s2[j-1])
            gap_left = matrix[i][j-1] +gapOpen
            gap_up = matrix[i-1][j]+gapOpen
            matrix[i][j] = max(0,match,gap_left,gap_up)
            if score <= matrix[i][j]:
                score = matrix[i][j]
                largest_position = (i,j)
    print("score {}, largest position: {}".format(score,largest_position))

    align = LocalAlign(matrix, s1, s2,gapOpen,largest_position[0],largest_position[1] )
    return matrix, align[0], align[1],score

def ALocal(s1,s2,gapOpen, gapExt):
    row = len(s1)+1
    col = len(s2)+1
    e = EMatrix(s1, s2,gapOpen,gapExt)
    f = FMatrix(s1, s2,gapOpen,gapExt)
    g = GMatrix(s1, s2,gapOpen,gapExt)
    v = VMatrix(s1, s2,gapOpen,gapExt)
    score = 0
    largest_position = (0,0)
    for i in range(1,row):
        for j in range(1,col):
            e[i][j] = max(e[i][j-1]+gapExt,v[i][j-1]+gapOpen+gapExt)
            f[i][j] = max(f[i-1][j]+gapExt,v[i-1][j]+gapOpen+gapExt)
            g[i][j] = v[i-1][j-1]+ Score(s1[i-1],s2[j-1])
            v[i][j] = max(0,g[i][j],e[i][j],f[i][j])
            if score <= v[i][j]:
                score = v[i][j]
                largest_position = (i,j)
    print("E:",e)
    print("F:",f)
    print("G:",g)
    print("V:",v)
    print(score)
    align = LocalAAlign(v, s1, s2,gapOpen, gapExt,largest_position[0],largest_position[1] )
    return v, align[0], align[1],score

#%%
def Main():
    print("Type the name of pattern file (should be located in same file with program code):")
    file = str(input())
    pattern = fileRead(file)
    while 1:
        print("Type 1 for Global Alignment with naive gap scoring \nType 2 for Global Alignment with affine gap scoring")
        print("Type 3 for Local Alignment with naive gap scoring \nType 4 for Local Alignment with affine gap scoring")
        print("For Exit type 0")
        mode = int(input())
        if mode == 0:
            break
        
        elif mode == 1:
            print("\nGap Opening Score:")
            gapOp = int(input())
            result = Global(pattern[1],pattern[0],gapOp)
            print(result)
            seq1 = result[1]
            seq2 = result[2]
            fileWrite(result[3],seq1,seq2)
            print("\nScore and Alignment is written to the result.txt\n")
            
        elif mode == 2:
            print("\nGap Opening Score:")
            gapOp = int(input())
            print("Gap Extension Score:")
            gapExt = int(input())
            result = AGlobal(pattern[1],pattern[0],gapOp,gapExt)
            print(result)
            seq1 = result[1]
            seq2 = result[2]
            fileWrite(result[3],seq1,seq2)
            print("\nScore and Alignment is written to the result.txt\n")
            
        elif mode == 3:
            print("\nGap Opening Score:")
            gapOp = int(input())
            result = Local(pattern[1],pattern[0],gapOp)
            print(result)
            seq1 = result[1]
            seq2 = result[2]
            fileWrite(result[3],seq1,seq2)
            print("\nScore and Alignment is written to the result.txt\n")
            
        elif mode == 4:
            print("\nGap Opening Score:")
            gapOp = int(input())
            print("Gap Extension Score:")
            gapExt = int(input())
            result = ALocal(pattern[1],pattern[0],gapOp,gapExt)
            print(result)
            seq1 = result[1]
            seq2 = result[2]
            fileWrite(result[3],seq1,seq2)
            print("\nScore and Alignment is written to the result.txt\n")
            
        else:
            print("You should enter a valid mode number\n")


Main()


