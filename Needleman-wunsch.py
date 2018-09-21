'''Implementasi Needleman-Wunsch Algorithm'''


import pandas as pd
import numpy as np

#Nilai Blosum62 diambil dari file csv.
blosum = pd.read_csv('blosum62.csv')

def needleman(seq1,seq2):
    def view_matrix(row_h, col_h, data, hide_zero=False, width=5):
        matrix = []
        cell = "%" + str(width) + "s"
        line = cell * (len(row_h) + 1)

        # print a header row
        matrix.append(line % tuple([' '] + list(row_h)))

        # print the data rows
        for b2, row in zip(col_h,data):
            if hide_zero:
                display = []
                for v in row:
                    if v == 0:
                        display.append('')
                    else:
                        display.append(v)
            else:
                display = row
            matrix.append(line % tuple([b2] + display))

        return '\n'.join(matrix)

    def score_matrix(seq1,seq2,blosum):
        score = []
        for a2 in seq2:
            current = []
            for a1 in seq1:
                current.append(blosum[a1][a2])
            score.append(current)
        return score

    score = score_matrix(seq1,seq2,blosum)
    print('\nMatrix Nilai Blosum')
    print(view_matrix(seq1,seq2,score))
    
    seq1 = " " + seq1
    seq2 = " " + seq2
    
    data = []
    for a in seq2:
        data.append(['-']*len(seq1))
    
    p=6
    data[0][0] = 0
    for i in range(1,len(seq2)):
        data[i][0] = data[i-1][0] - p

    for j in range(1, len(seq1)):
        data[0][j] = data[0][j-1] - p
    
    
    trace = []

    def total2(data,trace):
        gap = 6
        for p in seq2:
            trace.append(['-']*len(seq1))
        trace[0][0] = '*'
        
        for i in range(1, len(seq2)):
            for j in range(1, len(seq1)):
                diag = data[i-1][j-1] + score[i-1][j-1]
                ver = data[i-1][j] - gap
                hor = data[i][j-1] - gap
                data[i][j] = max(diag,ver,hor)
    
                if max(diag,hor,ver) == hor:
                    trace[i][j] = '-'
                elif max(diag,hor,ver) == ver:
                    trace[i][j] = '|'
                else:
                    trace[i][j] = '\ '
        return
    
    total2(data,trace)
    print('\nMatrix dengan Nilai Maksimal')
    print(view_matrix(seq1,seq2,data))
    print('\nMatrix Traceback')
    print(view_matrix(seq1,seq2,trace))
    
    align1 = ''
    align2 = ''
    
    
    def alignment(align1,align2,trace):
        current = None
        row = np.shape(trace)[0] - 1
        column = np.shape(trace)[1] - 1
    
        while current != '*':
            current = trace[row][column]
        
            if current == '\ ':
                align1 = seq1[column] + align1
                align2 = seq2[row] + align2
            
                row -= 1
                column -= 1
        
            elif current == '|':
                align1 = '-' + align1
                align2 = seq2[row] + align2
            
                row -= 1
        
            elif current == '-':
                align1 = seq1[column] + align1
                align2 = '-' + align2
            
                column -= 1
            elif current == '*':
                continue
            else:
                raise ValueError('Invalid value in traceback matrix: %s' % current)
            
            align_score=0
            for i in range(len(align1)):
                if (align1[i] == '-') | (align2[i] == '-') :
                    align_score += -6
                else:
                    align_score += blosum[align1[i]][align2[i]]
        return print('\n\nHasil Alignment: '),print(align1),print(align2), print('\nScore: ',align_score)
    
    
    
    return alignment(align1,align2,trace)

seq1 = input('Masukkan sequence protein pertama: ')
seq2 = input('Masukkan sequence protein kedua: ')

needleman(seq1,seq2)

