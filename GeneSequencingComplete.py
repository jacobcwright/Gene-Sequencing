#!/usr/bin/python3

from numpy import savez_compressed
from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT6':
	from PyQt6.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import math
import time
import random

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1

UP = 9
LEFT = 8
DIAG = 7

class GeneSequencing:

	def __init__( self ):
		pass
	
# This is the method called by the GUI.  _seq1_ and _seq2_ are two sequences to be aligned, _banded_ is a boolean that tells
# you whether you should compute a banded alignment or full alignment, and _align_length_ tells you 
# how many base pairs to use in computing the alignment

	def align( self, seq1, seq2, banded, align_length):
		# self.banded = banded
		# self.MaxCharactersToAlign = align_length
		if banded:
			score, alignment1, alignment2 = self.banded(seq1, seq2, align_length)
			return {'align_cost': score, 'seqi_first100': alignment1, 'seqj_first100': alignment2}
		else:
			score, alignment1, alignment2 = self.unrestricted(seq1, seq2, align_length)
			return {'align_cost': score, 'seqi_first100': alignment1, 'seqj_first100': alignment2}

	def unrestricted(self, seq1, seq2, align_length):
		# get dimensions of table
		row_end = min(len(seq1) + 1, align_length + 1)
		col_end = min(len(seq2) + 1, align_length + 1)

		# create dp table
		dp = []

		# iterate row by row and fill in table by getting min of 3 options, up, left, and diagonal
		for i in range(row_end):
			row = []
			for j in range(col_end):
				if i == 0:
					row.append((j * INDEL, LEFT))
				elif j == 0:
					row.append((i * INDEL, UP))
				else:
					diag_score = MATCH if seq1[i-1] == seq2[j-1] else SUB
					left = (row[j - 1][0] + INDEL, LEFT)
					up = (dp[i - 1][j][0] + INDEL, UP)
					diag = (dp[i - 1][j - 1][0] + diag_score, DIAG)
					row.append(min(left, up, diag, key=lambda x: x[0]))
			dp.append(row)

		score = dp[-1][-1][0]
		alignment1 = []
		alignment2 = []

		# backtrack to get alignment
		i = row_end - 1
		j = col_end - 1
		while i > 0 and j > 0:
			if dp[i][j][1] == UP:
				alignment1.append(seq1[i-1])
				alignment2.append('-')
				i -= 1
			elif dp[i][j][1] == LEFT:
				alignment1.append('-')
				alignment2.append(seq2[j-1])
				j -= 1
			elif dp[i][j][1] == DIAG:
				alignment1.append(seq1[i - 1])
				alignment2.append(seq2[j - 1])
				i -= 1
				j -= 1
		
		return score, "".join(alignment1[::-1]), "".join(alignment2[::-1])



	def banded(self, seq1, seq2, align_length):
		if abs(len(seq1) - len(seq2)) >= 500:
			return math.inf, "No Alignment Possible", "No Alignment Possible"
		banded_size = (2 * MAXINDELS) + 1
		band_end = min(len(seq1) + 1, align_length + 1)

		if len(seq1) > align_length:
			seq1 = seq1[:align_length]

		if len(seq2) > align_length:
			seq2 = seq2[:align_length]

		dp = []

		for i in range(band_end):
			row = [(math.inf, "")] * band_end
			start_col = i - MAXINDELS
			# end_col = min(len(seq2) + 1, align_length + 1, i + MAXINDELS + 1)
			for j in range(start_col, i + MAXINDELS + 1):
				if j < 0 or j > len(seq2) or j >= band_end:
					continue
				elif i == 0:
					row[j] = ((j * INDEL, LEFT))
					continue
				elif j == start_col and i < MAXINDELS:
					row[j] = ((i * INDEL, UP))
					continue
				else:
					if j == start_col:
						left = (math.inf, LEFT)
					else:
						left = (row[j - 1][0] + INDEL, LEFT)

					if j >= i + MAXINDELS:
						up = (math.inf, UP)
					else:
						up = (dp[i - 1][j][0] + INDEL, UP)

					diag_score = MATCH if seq1[i - 1] == seq2[j - 1] else SUB
					diag = (dp[i - 1][j - 1][0] + diag_score, DIAG)
					row[j] = (min(left, up, diag, key=lambda x: x[0]))
			dp.append(row)


		score = dp[-1][-1][0]
		alignment1 = []
		alignment2 = []

		# backtrack to get alignment
		i = len(dp) - 1
		j = len(dp[-1]) - 1
		while i > 0 and j > 0:
			if dp[i][j][1] == UP:
				alignment1.append(seq1[i-1])
				alignment2.append('-')
				i -= 1
			elif dp[i][j][1] == LEFT:
				alignment1.append('-')
				alignment2.append(seq2[j-1])
				j -= 1
			elif dp[i][j][1] == DIAG:
				alignment1.append(seq1[i - 1])
				alignment2.append(seq2[j - 1])
				i -= 1
				j -= 1
		
		return score, "".join(alignment1[::-1]), "".join(alignment2[::-1])
