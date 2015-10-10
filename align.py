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

# Computes global alignment and prints results
def runGlobalAlignment():
	alignments = findGlobalAligns()
	printResults(alignments)
	return

# Computes local alignment and prints results
def runLocalAlignment():
	alignments = findLocalAligns()
	printResults(alignments)
	return

# Find and prints out the highest global alignment score and
# returns a list of tuples with the coordinates of alignments with the highest score
def findGlobalAligns():
	# Stores the highest score
	max_score = 0
	# Stores coordinates of highest scoring alignments
	max_scorers = []
	# Scan last col of match matrix M
	for col in xrange(len(seq_a) + 1):
		curr_score = computeScores(len(seq_b), col, "M", False)
		if curr_score > max_score:
			# New max score
			max_scorers = [(len(seq_b), col)]
			max_score = curr_score
		elif round(curr_score,2) == round(max_score,2):
			# Add to max_scorers
			max_scorers.append((len(seq_b), col))
	# Scan last col of match matrix M
	for row in xrange(len(seq_b)):
		curr_score = computeScores(row, len(seq_a), "M", False)
		if curr_score > max_score:
			# New max score
			max_scorers = [(row, len(seq_a))]
			max_score = curr_score
		elif round(curr_score,2) == round(max_score,2):
			# Add to max_scorers
			max_scorers.append((row, len(seq_a)))
	# Write to file
	output_file.write(str(max_score) + "\n\n")
	print max_score
	print ""
	return max_scorers

# Finds and prints out the highest local alignment score and
# returns a list of tuples with the coordinates of alignments with the highest score
def findLocalAligns():
	# Stores the highest score
	max_score = 0
	# Stores coordinates of highest scoring alignments
	max_scorers = []
	# Scan entire match matrix M
	for col in xrange(len(seq_a) + 1):
		for row in xrange(len(seq_b)):
			curr_score = computeScores(row, col, "M", True)
			if curr_score > max_score:
				# New max score
				max_scorers = [(row, col)]
				max_score = curr_score
			elif round(curr_score,2) == round(max_score,2):
				# Add to max_scorers
				max_scorers.append((row, col))
	# Write to file
	output_file.write(str(max_score) + "\n\n")
	print max_score
	print ""
	return max_scorers

# Recursively computes and returns the alignment score of row i and col j in scoring matrix "curr_matrix"
# Populates all three scoring matrices with the scores and pointers up until row i, col j
# Can specify local or global alignments with the "local" argument
def computeScores(i, j, curr_matrix, local):
	# Base cases
	if i==0 or j==0:
		if curr_matrix == "M":
			M[i][j] = (0,[])
			return 0
		elif curr_matrix == "I_a":
			I_a[i][j] = (0,[])
			return 0
		else:
			I_b[i][j] = (0,[])
			return 0

	# Computes score for matrix M
	elif curr_matrix == "M":
		if type(M[i-1][j-1]) is not tuple:
			score = round(computeScores(i-1, j-1, "M", local) + S[alpha_b.index(seq_b[i-1])][alpha_a.index(seq_a[j-1])], 2)
			# Local alignment never stores negative scores
			if local is True and score < 0.0: m_score = 0
			else: m_score = score
		else:
			# Score was previously computed; saves computation
			score = round(M[i-1][j-1][0] + S[alpha_b.index(seq_b[i-1])][alpha_a.index(seq_a[j-1])], 2)
			# Local alignment never stores negative scores
			if local is True and score < 0.0: m_score = 0
			else: m_score = score

		if type(I_a[i-1][j-1]) is not tuple:
			score = round(computeScores(i-1, j-1, "I_a", local) + S[alpha_b.index(seq_b[i-1])][alpha_a.index(seq_a[j-1])], 2)
			if local is True and score < 0.0: ia_score = 0
			else: ia_score = score
		else:
			score = round(I_a[i-1][j-1][0] + S[alpha_b.index(seq_b[i-2])][alpha_a.index(seq_a[j-2])], 2)
			if local is True and score < 0.0: ia_score = 0
			else: ia_score = score

		if type(I_b[i-1][j-1]) is not tuple:
			score = round(computeScores(i-1, j-1, "I_b", local) + S[alpha_b.index(seq_b[i-1])][alpha_a.index(seq_a[j-1])], 2)
			if local is True and score < 0.0: ib_score = 0
			else: ib_score = score
		else:
			score = round(I_b[i-1][j-1][0] + S[alpha_b.index(seq_b[i-1])][alpha_a.index(seq_a[j-1])], 2)
			if local is True and score < 0.0: ib_score = 0
			else: ib_score = score

		# Find maximum score and create pointers to appropriate score matrices
		max_score = max(m_score, ia_score, ib_score)
		# If
		if max_score is 0 and M[i-1][j-1][1] == []:
			M[i][j] = (0, [])
		else:
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
		if type(M[i][j-1]) is not tuple:
			score = round(computeScores(i, j-1, "M", local) - d_b, 2)
			if local is True and score < 0.0: m_score = 0
			else: m_score = score
		else:
			score = round(M[i][j-1][0] - d_b, 2)
			if local is True and score < 0.0: m_score = 0
			else: m_score = score

		if type(I_a[i][j-1]) is not tuple:
			score = round(computeScores(i, j-1, "I_a", local) - e_b, 2)
			if local is True and score < 0.0: ia_score = 0
			else: ia_score = score
		else:
			score = round(I_a[i][j-1][0] - e_b, 2)
			if local is True and score < 0.0: ia_score = 0
			else: ia_score = score

		# Find maximum score and create pointers to appropriate score matrices
		max_score = max(m_score, ia_score)
		if max_score is 0 and M[i-1][j-1][1] == []:
			I_a[i][j] = (0, [])
		else:
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
		if type(M[i-1][j]) is not tuple:
			score = round(computeScores(i-1, j, "M", local) - d_a, 2)
			if local is True and score < 0.0: m_score = 0
			else: m_score = score
		else:
			score = round(M[i-1][j][0] - d_a, 2)
			if local is True and score < 0.0: m_score = 0
			else: m_score = score

		if type(I_b[i-1][j]) is not tuple:
			score = round(computeScores(i-1, j, "I_b", local) - e_a, 2)
			if local is True and score < 0.0: ib_score = 0
			else: ib_score = score
		else:
			score = round(I_b[i-1][j][0] - e_a, 2)
			if local is True and score < 0.0: ib_score = 0
			else: ib_score = score

		# Find maximum score and create pointers to appropriate score matrices
		max_score = max(m_score, ib_score)
		if max_score is 0 and M[i-1][j-1][1] == []:
			I_b[i][j] = (0, [])
		else:
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

def printResults(alignments):
	for startpoint in alignments:
		printAlignment(startpoint, "M", "", "")
	f.close()
	output_file.close()

def printAlignment(point, matrix, align_a, align_b):
	row = point[0]
	col = point[1]
	if row is 0 or col is 0:
		output_file.writelines(str(align_a)+'\n'+str(align_b)+'\n\n')
		print align_a
		print align_b
		print ""
	elif matrix == "M":
		if M[row][col][1] == []:
			output_file.writelines(str(align_a)+'\n'+str(align_b)+'\n\n')
			print align_a
			print align_b
			print ""
		else:
			for path in M[row][col][1]:
				printAlignment( (row-1, col-1), path, seq_a[col-1] + align_a, seq_b[row-1] + align_b)
	elif matrix == "I_a":
		if M[row][col][1] == []:
			output_file.writelines(str(align_a)+'\n'+str(align_b)+'\n\n')
			print align_a
			print align_b
			print ""
		else:
			for path in I_a[row][col][1]:
				printAlignment( (row, col-1), path, seq_a[col-1] + align_a, "_" + align_b)
	elif matrix == "I_b":
		if M[row][col][1] == []:
			output_file.writelines(str(align_a)+'\n'+str(align_b)+'\n\n')
			print align_a
			print align_b
			print ""
		else:
			for path in I_b[row][col][1]:
				printAlignment( (row-1, col), path, "_" + align_a, seq_b[row-1] + align_b)

if int(input_file[2].strip()) == 0:
	runGlobalAlignment()
	print "Global"
else:
	runLocalAlignment()
	print "Local"

"""
for row in xrange(8):
	string=[]
	for col in xrange(5):
		string += I_b[97+row][col]
	print string
	"""