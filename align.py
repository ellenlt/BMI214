# Aligns two DNA sequences with affine gap penalty, either using global ends-free or local alignment.
# Accepts two arguments, the input_file and output_file
# i = row (sequence B), j = col (sequence A)

import sys
print "Input file: '%s'" % sys.argv[1]
print "Output file: '%s'" % sys.argv[2]

# Read in input_file
# Store sequences A and B as seq_a and seq_B
# Initialize empty score matrices M, Ia, Ib where width = seq_a length + 1 and v.v.
# Determine local (1) or global (0) alignment
# Stores gap penalties
# Stores alphabet order in alpha_a and alpha_b
# Create match matrix S(a,b) where a = col / index along alphabet A and b = row / index along alphabet B

f = open(sys.argv[1])
input_file = f.readlines()

# Original sequences for alignment
seq_a = input_file[0].strip()
seq_b = input_file[1].strip()

#if input_file[2] == 0:
	# run global alignment
	# else:
	# run local alignment

penalties = input_file[3].split(" ")
d_a = int(penalties[0])	# Gap open penalty for sequence A
e_a = int(penalties[1])	# Gap extension penalty for sequence A
d_b = int(penalties[2])	# Gap open penalty for sequence B
e_b = int(penalties[3])	# Gap extension penalty for sequence B

alpha_a_size = int(input_file[4])
alpha_a = input_file[5].strip()	# Alphabet used for sequence A
alpha_b_size = int(input_file[6])
alpha_b = input_file[7].strip()	# Alphabet used for sequence B

# Initialize scoring matrices M, I_a, and I_b to initial values of 0
# Sequence A along top (spans columns), Sequence B along left (spans rows)
# Each cell contains an alignment score and may also contain string pointers to the matrix/matrices that produced the score
M = [[0 for a in xrange(len(seq_a)+1)] for b in xrange(len(seq_b)+1)]
I_a = [[0 for a in xrange(len(seq_a)+1)] for b in xrange(len(seq_b)+1)]
I_b = [[0 for a in xrange(len(seq_a)+1)] for b in xrange(len(seq_b)+1)]

# Set first row and column to 0
"""
for i in xrange(len(seq_a)+1):
	M[0][i] = 0
	I_a[0][i] = 0
	I_b[0][i] = 0
for j in xrange(len(seq_b)+1):
	M[j][0]= 0
	I_a[j][0] = 0
	I_b[j][0] = 0
"""
# Initialize and populate match matrix S
# Alphabet A along top (spans columns), Alphabet B along left (spans rows)
# Tells you how many points to assign to a match/mismatch
S = [[0 for a in xrange(alpha_a_size)] for b in xrange(alpha_b_size)]
for i in xrange(alpha_a_size*alpha_b_size):
	entry = input_file[8+i].split(" ")
	S[int(entry[1])-1][int(entry[0])-1] = int(entry[4])

def runGlobalAlignment():
	alignments = findGlobalAligns()
#	printResults(alignments)
	return

def findGlobalAligns():
	computeScores(len(seq_b), len(seq_a), "M")
	max_score = 0
	max_scorers = []
	# Scan last row for highest score
	for col in xrange(alpha_a_size + 1)				#implement this so that you can compute all last col last row and figure out
		curr_score = computeScore(i, alpha_b_size, "M")  #highest scorers. then, check if error cell is correct
		if curr_score > max_score:	# New maximum score
			max_scorers = [(i, alpha_b_size)]
			max_score = curr_score
		elif curr_score == max_score:	# Add to maximum scorer
			max_scorers.append(i, alpha_b_size)
	# Scan last column for highest score, not including bottom right which was already scanned
	for j in xrange(alpha_b_size + 1)
		curr_score = computeScore(alpha_a_size, j, "M")
		if curr_score > max_score:
			# replace max_scorers
			max_scorers = [(alpha_a_size, j)]
			max_score = curr_score
		elif curr_score == max_score:
			# add to max_scorers
			max_scorers.append(alpha_a_size, j)
	return max_scorers

# Recursively computes global alignment scoring matrices up until row i, col j
def computeScores(i, j, curr_matrix):
	# Base cases
#	print "at first, (i,j) = (%d,%d) and m = %s" % (i, j, curr_matrix)
	if i==0 or j==0:
#		print "we set something -> zero for m = %s" % curr_matrix
		if curr_matrix == "M":
			M[i][j] = 0
			return 0
		elif curr_matrix == "I_a":
			I_a[i][j] = 0
			return 0
		else:
			I_b[i][j] = 0
			return 0
	elif curr_matrix == "M":
		# Compute scores for M, I_a, and I_b
		if type(M[i-1][j-1]) is not tuple: m_score = computeScores(i-1, j-1, "M") + S[alpha_b.index(seq_b[i-1])][alpha_a.index(seq_a[j-1])]
		else: m_score = M[i-1][j-1][0] + S[alpha_b.index(seq_b[i-1])][alpha_a.index(seq_a[j-1])]
		print S[alpha_b.index(seq_b[i-1])][alpha_a.index(seq_a[j-1])]

		if type(I_a[i-1][j-1]) is not tuple: ia_score = computeScores(i-1, j-1, "I_a") + S[alpha_b.index(seq_b[i-1])][alpha_a.index(seq_a[j-1])]
		else: ia_score = I_a[i-1][j-1][0] + S[alpha_b.index(seq_b[i-2])][alpha_a.index(seq_a[j-2])]

		if type(I_b[i-1][j-1]) is not tuple: ib_score = computeScores(i-1, j-1, "I_b") + S[alpha_b.index(seq_b[i-1])][alpha_a.index(seq_a[j-1])]
		else: ib_score = I_b[i-1][j-1][0] + S[alpha_b.index(seq_b[i-1])][alpha_a.index(seq_a[j-1])]

		# Find maximum score and create pointers to appropriate score matrices
		max_score = max(m_score, ia_score, ib_score)
		max_scorers = []
		if(m_score == max_score): max_scorers.append("M")
		if(ia_score == max_score): max_scorers.append("I_a")
		if(ib_score == max_score): max_scorers.append("I_b")

		# Store the score and pointers
#		print "finally, (i,j) = (%d,%d) and m = %s" % (i, j, curr_matrix)
		M[i][j] = (max_score, max_scorers)
		return max_score

	elif curr_matrix == "I_a":
		# Compute scores for M, I_a, and I_b
		if type(M[i][j-1]) is not tuple: m_score = computeScores(i, j-1, "M") - d_b
		else: m_score = M[i][j-1][0] - d_b

		if type(I_a[i][j-1]) is not tuple: ia_score = computeScores(i, j-1, "I_a") - e_b
		else: ia_score = I_a[i][j-1][0] - e_b

		# Find maximum score and create pointers to appropriate score matrices
		max_score = max(m_score, ia_score)
		max_scorers = []
		if(m_score == max_score): max_scorers.append("M")
		if(ia_score == max_score): max_scorers.append("I_a")

		# Store the score and pointers
		I_a[i][j] = (max_score, max_scorers)
		return max_score
		
	else:
		# Compute scores for M, I_a, and I_b
		if type(M[i-1][j]) is not tuple: m_score = computeScores(i-1, j, "M") - d_a
		else: m_score = M[i-1][j][0] - d_a

		if type(I_b[i-1][j]) is not tuple: ib_score = computeScores(i-1, j, "I_b") - e_a
		else: ib_score = I_b[i-1][j][0] - e_a

		# Find maximum score and create pointers to appropriate score matrices
		max_score = max(m_score, ib_score)
		max_scorers = []
		if(m_score == max_score): max_scorers.append("M")
		if(ib_score == max_score): max_scorers.append("I_b")

		# Store the score and pointers
		I_b[i][j] = (max_score, max_scorers)
		return max_score

	return

# trace_from_gap_a(i, j)
# i = # index of current letter in seq_a, j = # index of current letter in seq_b
# Base case: If i=0 or j=0, return 0
# Otherwise, return the maximum of the following:
# trace_from_match(i-1, j) + d_y
# trace_from_gap_a(i-1, j) - e_y

# trace_from_gap_b(i, j)
# i = # index of current letter in seq_a, j = # index of current letter in seq_b
# Base case: If i=0 or j=0, return 0
# Otherwise, return the maximum of the following:
# trace_from_match(i, j-1) + d_x
# trace_from_gap_b(i, j-1) - e_x

def runLocalAlignment():
	# scan entire matrix to find best score
	# Never record a negative score in score matrix
	return

def printResults(alignments):
	aligned_a = ""
	aligned_b = ""
	# Use the row:col tuples in alignments and the M matrix to traceback
	# Construct string
	return

output_file = open(sys.argv[2], 'w')
# Write file to output_file
# Print file to output


runGlobalAlignment()
print "M:"
for row in range(len(M)):
	print M[row]
print "Ia:"
for row in range(len(M)):
	print I_a[row]
print "Ib:"
for row in range(len(M)):
	print I_b[row]