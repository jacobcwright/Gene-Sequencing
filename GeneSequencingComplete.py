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

class GeneSequencing:
	score = 0
	alignment1 = ''
	alignment2 = ''

	def __init__( self ):
		pass
	
# This is the method called by the GUI.  _seq1_ and _seq2_ are two sequences to be aligned, _banded_ is a boolean that tells
# you whether you should compute a banded alignment or full alignment, and _align_length_ tells you 
# how many base pairs to use in computing the alignment

	def align( self, seq1, seq2, banded, align_length):
		self.banded = banded
		self.MaxCharactersToAlign = align_length
		return {'align_cost':score, 'seqi_first100':alignment1, 'seqj_first100':alignment2}

	def unrestricted(self, seq1, seq2, align_length):
		# create dp table
		dp = []
		# iterate row by row and fill in table by getting min of 3 options, up, left, and diagonal
		# add a tuple into the table index with the value and the direction it came from

		# then add row to table

		# the last index of dp is the answer (score)
		# then trace back to get the alignment (alignment1 and alignment2)



	def banded(self, seq1, seq2, align_length):
		return