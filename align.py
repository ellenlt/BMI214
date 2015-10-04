import sys

# Aligns two DNA sequences with affine gap penalty, either using global ends-free or local alignment.
# Accepts two arguments, the input_file and output_file
# i = row (sequence B), j = col (sequence A)
print "Input file: '%s'" % sys.argv[1]
print "Output file: '%s'" % sys.argv[2]

# Read in input_file
f = open(sys.argv[1])
input_file = f.readlines()

# Create new output file
output_file = open(sys.argv[2], 'w')

# Original sequences for alignment
seq_a = input_file[0].strip()
seq_b = input_file[1].strip()

#if input_file[2] == 0:
	# run global alignment
	# else:
	# run local alignment

penalties = input_file[3].split(" ")
d_a = float(penalties[0])	# Gap open penalty for sequence A
e_a = float(penalties[1])	# Gap extension penalty for sequence A
d_b = float(penalties[2])	# Gap open penalty for sequence B
e_b = float(penalties[3])	# Gap extension penalty for sequence B

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

# Initialize and populate match matrix S
# Alphabet A along top (spans columns), Alphabet B along left (spans rows)
# Tells you how many points to assign to a match/mismatch
S = [[0 for a in xrange(int(alpha_a_size))] for b in xrange(int(alpha_b_size))]
for i in xrange(int(alpha_a_size)*int(alpha_b_size)):
	entry = input_file[8+i].split()
	S[int(entry[1])-1][int(entry[0])-1] = (float(entry[4]))

def runGlobalAlignment():
	alignments = findGlobalAligns()
	printResults(alignments)
	return

# Scans the last row and column to find the highest scores.
# Prints out the highest score and
# returns a list of tuples with the coordinates of the highest scores
def findGlobalAligns():
	max_score = 0
	max_scorers = []
	# Scan last row and col for highest score and stores the coordinates in max_scorers
	for col in xrange(len(seq_a) + 1):
		curr_score = computeScores(len(seq_b), col, "M")
		if curr_score > max_score:				# New maximum score
			max_scorers = [(len(seq_b), col)]
			max_score = curr_score
		elif round(curr_score,2) == round(max_score,2):			# Add to maximum scorer
			max_scorers.append((len(seq_b), col))
	# Scan last column for highest score, not including bottom right which was already scanned
	for row in xrange(len(seq_b)):
		curr_score = computeScores(row, len(seq_a), "M")
		if curr_score > max_score:
			# replace max_scorers
			max_scorers = [(row, len(seq_a))]
			max_score = curr_score
		elif round(curr_score,2) == round(max_score,2):
			# add to max_scorers
			max_scorers.append((row, len(seq_a)))
	# Write to file
	output_file.write(str(max_score) + "\n\n")
	print max_score
	print ""
	return max_scorers

# Returns the recursively computed global alignment score of row i, col j
# and populates all three scoring matrices with the scores and pointers up until row i, col j
def computeScores(i, j, curr_matrix):
	# Base cases
	if i==0 or j==0:
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
		if type(M[i-1][j-1]) is not tuple: m_score = round(computeScores(i-1, j-1, "M") + S[alpha_b.index(seq_b[i-1])][alpha_a.index(seq_a[j-1])], 2)
		else:
			m_score = round(M[i-1][j-1][0] + S[alpha_b.index(seq_b[i-1])][alpha_a.index(seq_a[j-1])], 2)

		if type(I_a[i-1][j-1]) is not tuple: ia_score = round(computeScores(i-1, j-1, "I_a") + S[alpha_b.index(seq_b[i-1])][alpha_a.index(seq_a[j-1])], 2)
		else: ia_score = round(I_a[i-1][j-1][0] + S[alpha_b.index(seq_b[i-2])][alpha_a.index(seq_a[j-2])], 2)

		if type(I_b[i-1][j-1]) is not tuple: ib_score = round(computeScores(i-1, j-1, "I_b") + S[alpha_b.index(seq_b[i-1])][alpha_a.index(seq_a[j-1])], 2)
		else: ib_score = round(I_b[i-1][j-1][0] + S[alpha_b.index(seq_b[i-1])][alpha_a.index(seq_a[j-1])], 2)

		# Find maximum score and create pointers to appropriate score matrices
		max_score = max(m_score, ia_score, ib_score)
		max_scorers = []
		if(round(m_score,2) == round(max_score,2)): max_scorers.append("M")
		# If pointers trace back to a 0 in both the M and the I_a or I_b matrices, only store the pointer to the M matrix
		if(round(ia_score,2) == round(max_score,2)):
			if(type(I_a[i-1][j-1]) is tuple and I_a[i-1][j-1][0] != 0):
				max_scorers.append("I_a")
		if(round(ib_score,2) == round(max_score,2) and I_b[i-1][j-1] != 0):
			if(type(I_b[i-1][j-1]) is tuple and I_b[i-1][j-1][0] != 0):
				max_scorers.append("I_b")

		# Store the score and pointers
		M[i][j] = (max_score, max_scorers)
		return max_score

	elif curr_matrix == "I_a":
		# Compute scores for M, I_a, and I_b
		if type(M[i][j-1]) is not tuple: m_score = round(computeScores(i, j-1, "M") - d_b, 2)
		else: m_score = round(M[i][j-1][0] - d_b, 2)

		if type(I_a[i][j-1]) is not tuple: ia_score = round(computeScores(i, j-1, "I_a") - e_b, 2)
		else: ia_score = round(I_a[i][j-1][0] - e_b, 2)

		# Find maximum score and create pointers to appropriate score matrices
		max_score = max(m_score, ia_score)
		max_scorers = []
		if(round(m_score,2) == round(max_score,2)): max_scorers.append("M")
		# If pointers trace back to a 0 in both the M and the I_a or I_b matrices, only store the pointer to the M matrix
		if(round(ia_score,2) == round(max_score,2)):
			if(type(I_a[i][j-1]) is tuple and I_a[i][j-1][0] != 0): max_scorers.append("I_a")

		# Store the score and pointers
		I_a[i][j] = (max_score, max_scorers)
		return max_score
		
	else:
		# Compute scores for M, I_a, and I_b
		if type(M[i-1][j]) is not tuple: m_score = round(computeScores(i-1, j, "M") - d_a, 2)
		else: m_score = round(M[i-1][j][0] - d_a, 2)

		if type(I_b[i-1][j]) is not tuple: ib_score = round(computeScores(i-1, j, "I_b") - e_a, 2)
		else: ib_score = round(I_b[i-1][j][0] - e_a, 2)

		# Find maximum score and create pointers to appropriate score matrices
		max_score = max(m_score, ib_score)
		max_scorers = []
		if(round(m_score,2) == round(max_score,2)):
			max_scorers.append("M")
		# If pointers trace back to a 0 in both the M and the I_a or I_b matrices, only store the pointer to the M matrix
		if(round(ib_score,2) == round(max_score,2)): 
			if(type(I_b[i-1][j]) is tuple and I_b[i-1][j][0] != 0):
				max_scorers.append("I_b")

		# Store the score and pointers
		I_b[i][j] = (max_score, max_scorers)
		return max_score

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

#For first letter, check if 

def printResults(alignments):
	for startpoint in alignments:
		printAlignment(startpoint, "M", "", "")

def printAlignment(point, matrix, align_a, align_b):
	row = point[0]
	col = point[1]
	if row is 0 or col is 0:
		output_file.writelines(str(align_a)+'\n'+str(align_b)+'\n\n')
		print align_a
		print align_b
		print ""
	elif matrix == "M":
		for path in M[row][col][1]:
			printAlignment( (row-1, col-1), path, seq_a[col-1] + align_a, seq_b[row-1] + align_b)
	elif matrix == "I_a":
		for path in I_a[row][col][1]:
			printAlignment( (row, col-1), path, seq_a[col-1] + align_a, "_" + align_b)
	elif matrix == "I_b":
		for path in I_b[row][col][1]:
			printAlignment( (row-1, col), path, "_" + align_a, seq_b[row-1] + align_b)
	
runGlobalAlignment()