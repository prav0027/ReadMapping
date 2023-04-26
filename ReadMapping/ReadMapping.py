class BWA:
	#this initializer function creates all the datastructures necessary using the reference string
	def __init__(self, reference):
		#declare datastructures
		rotation_list, rotation_list_reverse, suffix_array, bwt = [list() for i in range(4)]
		C, Occ, Occ_reverse = [dict() for i in range(3)]
		alphabet = set()
		reverse_reference = reference[::-1]#reverse reference
		
		#Construct the alphabet. (This would be hard coded for DNA examples)
		reference = reference.lower()
		for char in reference:
			alphabet.add(char)
		
		#initialize 2 auxillary datastructures
		for char in alphabet:
			C[char] = 0
			Occ[char] = list()# in Occ, each character has an associated list of integer values (for each index along the reference)
			Occ_reverse[char] = list()
	
		#append the ending character to the reference string
		reference = "%s$" % reference
		reverse_reference = "%s$" % reverse_reference

		#create all the rotation/suffix combinations of the reference and reverse reference, and their starting index positions
		for i in range(len(reference)):
			new_rotation = "%s%s" % (reference[i:],reference[0:i])
			struct = Suffix(new_rotation,i)
			rotation_list.append(struct)
			
			new_rotation_reverse = "%s%s" % (reverse_reference[i:],reverse_reference[0:i])
			struct_rev = Suffix(new_rotation_reverse,i)
			rotation_list_reverse.append(struct_rev)
		
			#create the C datastructure. C(a) = the number of characters 'a' in the Reference that are lexographically smaller than 'a'
			#NOTE, the C datastructure is not required for the reverse reference
			if reference[i]!='$':
				for char in alphabet:
					if reference[i] < char:
						C[char] = C[char] + 1	
		
		#sort the rotations/suffixes using the suffix/rotation text as the key
		rotation_list.sort(key=textKey)
		rotation_list_reverse.sort(key=textKey)
	
		#now record the results into 2 seperate lists, the suffix (or S) array and the BWT (or B) array
		#also calculate the auxilliary datastructure Occ (or O)
		for i in rotation_list:
			suffix_array.append(i.pos)#the position of the reordered suffixes forms the Suffix Array elements
			bwt.append(i.text[-1:])#the last character in each rotation (in the new order) forms the BWT string elements
		
			#now construct the Occ (or C) datastructure
			for char in alphabet:
				if len(Occ[char]) == 0:
					prev = 0
				else:
					prev = Occ[char][-1]
				if i.text[-1:] == char:
					Occ[char].append(prev+1)
				else:
					Occ[char].append(prev)
					
		#now record the results into 2 seperate lists, the suffix (or S) array and the BWT (or B) array
		#also calculate the auxilliary datastructures, C and Occ (or O)
		for i in rotation_list_reverse:
			#construct the Occ (or C) datastructure
			for char in alphabet:
				if len(Occ_reverse[char]) == 0:
					prev = 0
				else:
					prev = Occ_reverse[char][-1]
				if i.text[-1:] == char:
					Occ_reverse[char].append(prev+1)
				else:
					Occ_reverse[char].append(prev)					
					
		#save all the useful datastructures as class variables for easy future access
		self.SA = suffix_array
		self.BWT = bwt
		self.C = C
		self.Occ = Occ
		self.Occ_reverse = Occ_reverse #the Occ datastructure for the reverse reference, using to construct the D array (the lower bound on the number of differences allowed), to speed up alignments 
		self.n = len(reference)
		self.D = list()#empty list for later use
		self.alphabet = alphabet

	#get the position(s) of the query in the reference
	def find_match(self,query,num_differences):
		if num_differences == 0:
			return self.exact_match(query)
		else:
			return self.inexact_match(query,num_differences)

	#exact matching - no indels or mismatches allowed
	def exact_match(self, query):
		query = query.lower()
		i = 0
		j = self.n - 1
		
		for x in range(len(query)):
			newChar = query[-x-1]
			newI = self.C[newChar] + self.OCC(newChar,i-1) + 1
			newJ = self.C[newChar] + self.OCC(newChar,j)
			i = newI
			j = newJ
		matches = self.SA[i:j+1]
		return matches

	#inexact matching, z is the max threshold for allowed edits
	def inexact_match(self,query,z):
		self.calculate_d(query)
		SA_indeces = self.inexact_recursion(query, len(query)-1, z, 0, self.n-1)
		return [self.SA[x] for x in SA_indeces]#return the values in the SA

	#recursion function that effectively "walks" through the suffix tree using the SA, BWT, Occ and C datastructures
	def inexact_recursion(self,query,i,z,k,l):
		tempset = set()
		resultZ = dict()
		resultDict = dict()
		insertionList = list()
		#2 stop conditions, one when too many differences have been encountered, another when the entire query has been matched, terminating in success
		#if (z < self.get_D(i) and use_lower_bound_tree_pruning) or (z < 0 and not use_lower_bound_tree_pruning):#reached the limit of differences at this stage, terminate this path traversal
		if (z < self.get_D(i)):
			#if debug:print ("too many differences, terminating path\n" )
			return set()
			#return set(), insertionList #return empty set	
		if i < 0:#empty query string, entire query has been matched, return SA indexes k:l
			#if debug:print ("query string finished, terminating path, success! k=%d, l=%d\n" % (k,l))
			for m in range(k,l+1):
				#print("k: ", k, "; l: ", l, "; m: ", m)
				tempset.add(m)
				resultZ[m] = z
				resultDict[m] = resultZ
			return tempset#, insertionList
			return resultDict
			
		result = set()
		#if indels_allowed:
		resultAdd = self.inexact_recursion(query,i-1,z-insertion_penalty,k,l)
		result = result.union(resultAdd) #without finding a match or altering k or l, move on down the query string. Insertion
			#if len(insertionAdd)>0:
			#	print("insertionAdd insert: ", insertionAdd)
			#	insertionList.append(insertionAdd)
			#insertionList.append(i)
		for char in self.alphabet:#for each character in the alphabet
			#find the SA interval for the char
			newK = self.C[char] + self.OCC(char,k-1) + 1 
			newL = self.C[char] + self.OCC(char,l)
			if newK <= newL:#if the substring was found
				#if indels_allowed:
					#print("deletion")
				resultAdd = self.inexact_recursion(query,i,z-deletion_penalty,newK,newL)
				result = result.union(resultAdd)
					#if len(insertionAdd)>0:
					#	print("insertionAdd deletion: ", insertionAdd)
					#	insertionList.append(insertionAdd)
					#result = result.union(self.inexact_recursion(query,i,z-deletion_penalty,newK,newL))# Deletion
				#if debug:print ("char '%s found' with k=%d, l=%d. z = %d: parent k=%d, l=%d" % (char,newK,newL,z,k,l))
				if char == query[i]:#if the char was correctly aligned, then continue without decrementing z (differences)
					#print("char is correctly aligned")
					#resultAdd, insertionAdd = self.inexact_recursion(query,i-1,z,newK,newL)
					#result = result.union(resultAdd)
					#if len(insertionAdd)>0:
					#	print("insertionAdd correctly aligned: ", insertionAdd)
					#	insertionList.append(insertionAdd)
					#	print("insertionList: ", insertionList)
					result = result.union(self.inexact_recursion(query,i-1,z,newK,newL))
				else:#continue but decrement z, to indicate that this was a difference/unalignment
					#print("there is a difference")
					#resultAdd, insertionAdd = self.inexact_recursion(query,i-1,z-mismatch_penalty,newK,newL)
					#result = result.union(resultAdd)
					#if len(insertionAdd)>0:
					#	print("insertionAdd difference: ", insertionAdd)
					#	insertionList.append(insertionAdd)
					result = result.union(self.inexact_recursion(query,i-1,z-mismatch_penalty,newK,newL))
		#print("result: ", result, 'z: ', z)
		return result#, insertionList

	#calculates the D array for a query, used to prune the tree walk and increase speed for inexact searching
	def calculate_d(self,query):
		k = 0
		l = self.n-1
		z = 0
		self.D = list()#empty the D array
		for i in range(len(query)):
			k = self.C[query[i]] + self.OCC(query[i],k-1,reverse=True) + 1
			l = self.C[query[i]] + self.OCC(query[i],l,reverse=True)
			if k > l:#if this character has NOT been found
				k = 0
				l = self.n - 1
				z = z + 1
			self.D.append(z)

	#returns normal Occ value, otherwise returns the reverse Occ if explicitly passed as an argument
	#NOTE Occ('a',-1) = 0 for all 'a'
	def OCC(self,char,index,reverse=False):
		if index < 0:
			return 0
		else:
			if reverse:
				return self.Occ_reverse[char][index]
			else:
				return self.Occ[char][index]
	
	#gets values from the D array
	#NOTE D(-1) = 0
	def get_D(self,index):
		if index < 0:
			return 0
		else:
			return self.D[index]

class Suffix:
	"""
	A simple class with 2 variables, used for sorting and calculating the Suffix Array and BWT array
	Each instance holds the position of the suffix and the suffix (text) itself
	"""
	def __init__(self, text, position):
		self.text = text
		self.pos = position

#this is used to sort the Suffix objects, according to their text key
def textKey( a ): return a.text

insertionsDict = dict()
deletionsDict = dict()
substitutionsDict = dict()

def findMutations(query, index, numberOfMismatches, mismatchesCount, data):
	#mismatchesCount=0
	for i in range(len(query)):
		substitutionNeeded = True
		#print("i = ", i, ", query[i]: ", query[i])
		#print("index = ", index, ", reference[index]: ", reference[index])
		#print("before insertion numberofmismatches: ", numberOfMismatches)
		mismatchesCountTemp = mismatchesCount
		insertions_i_list = []
		insertions_char_list = []
		deletions_i_list = []
		if reference[index]==query[i]:
			index+=1
			continue
		#print("i = ", i, ", query[i]: ", query[i])
		#print("index = ", index, ", reference[index]: ", reference[index])
		#check for insertions
		while numberOfMismatches>0:#mismatchesCountTemp < numberOfMismatches and mismatchesCount < numberOfMismatches:
			iTemp = i+1
			if index in data.find_match(query[iTemp:], numberOfMismatches-1):
				print("insertion @ index: ", index)
				substitutionNeeded = False
				insertions_char_list.append(query[i])
				insertions_i_list.append(index-1) #switch from i to correctindex
				index+=-1
				numberOfMismatches+=-1
				mismatchesCountTemp+=1
				mismatchesCount = mismatchesCountTemp
			else:
				#if len(insertions_i_list)>0:
				#	pop = insertions_char_list.pop()
				#	print("get rid of insertion @ index: ", pop)
				#	insertions_i_list.pop()
				#	numberOfMismatches+=1
				#	mismatchesCount+=-1
				break
		for i in range(len(insertions_i_list)):
			if insertions_i_list[i] not in list(insertionsDict.keys()):
				insertionsDict[insertions_i_list[i]] = list(insertions_char_list[i])
			else:
				#print(insertionsDict[insertions_i_list[i]])
				tempList = insertionsDict[insertions_i_list[i]]
				tempList.append(insertions_char_list[i])
				insertionsDict[insertions_i_list[i]] = tempList
			#if i is a key in insertionDict, then add each i to the existing list (value)
			#else, add 
			#insertionsDict[insertions_i_list[i]] = insertions_char_list[i]
		#print("insertionsDict: ", insertionsDict)
		#print("number of mismatches: ", numberOfMismatches)
		#print("after insertion numberofmismatches: ", numberOfMismatches)
		while numberOfMismatches>0:#mismatchesCountTemp < numberOfMismatches and mismatchesCount < numberOfMismatches:
			if index+1 in data.find_match(query[i:], numberOfMismatches-1):
				#print("can do deletion at i: ", i)
				substitutionNeeded = False
				deletions_i_list.append(index)
				#print("deletions: ", index)
				index+=1
				numberOfMismatches+=-1
				mismatchesCountTemp+=1
				mismatchesCount = mismatchesCountTemp
			else:
				#print("cannot do deletion, ", query[i:], " i: ", i, " number of mismatches - 1: ", numberOfMismatches-1)
				#if len(deletions_i_list)>0:
				#	pop = deletions_i_list.pop()
				#	print("pop: ", pop)
				#	numberOfMismatches+=1
				#	mismatchesCount+=-1
				break
		#print("deletions list: ", deletions_i_list)
		for val in deletions_i_list:
			if val not in list(deletionsDict.keys()):
				deletionsDict[val] = list(reference[val])
			else:
				#print(val)
				#print(deletionsDict)
				tempList = deletionsDict[val]
				tempList.append(reference[val])
				deletionsDict[val] = tempList
			#deletionsDict[val] = reference[val]

		#print("deletionsDict: ", deletionsDict)
		#print("after deletion numberofmismatches: ", numberOfMismatches)
		if substitutionNeeded==True:
			#index_i_tuple = (index, i)
			nucleotide_tuple = (reference[index], query[i])
			numberOfMismatches+=-1
			#print("performing substitution at index: ", index, " and i: ", i)
			if index not in list(substitutionsDict.keys()):
				substitutionTempList = []
				substitutionTempList.append(nucleotide_tuple)
				substitutionsDict[index] = substitutionTempList
			else:
				tempList = substitutionsDict[index]
				tempList.append(nucleotide_tuple)
				substitutionsDict[index] = tempList
		index+=1
		#print("after substitution numberofmismatches: ", numberOfMismatches)
	return 

def read_reads(read_fn):
    all_reads = []
    with open(read_fn, 'r') as f:
        #next(f)
        for line in f:
            if '>' in line:
                continue
            line = line.strip()
            all_reads.append(line)
    return all_reads

def read_reference(ref_fn):
    with open(ref_fn, 'r') as f:
        next(f)
        output_reference = ''
        for line in f:
            if '>' in line:
                continue
            line = line.strip()
            output_reference += line  # We append each line to the output reference string.
    return output_reference

#environment variables
#debug = True
#use_lower_bound_tree_pruning = True #set this to false (in conjunction with debug=True) to see the full search through the suffix trie
#search parameters
#indels_allowed = True # turn off for mismatches only, no insertion or deletions allowed
insertion_penalty = 1
deletion_penalty = 1
mismatch_penalty = 1

#reference and query strings
#reference = """TTATAAACACAGGACGAGCGCTCCGGATCAAAAACAACCAGTCTGGCTAAACGAGTAACTCGACCCCGAGTGTGAGCAATCGTAGACGTCTGTGGTATTGGGCAAAGGTTTTAGAAATTGCTATGGGCCCTATAGTCATTTGGGGCTTGCTCCTATAGTTCTCCGTATCCAGTTGTGCTAATGGGAGGTCGCCAGGCGGGGGACCAACTATGCCCCACAGGACAAATCTGACGCCGTGATTGCAGCCCACAAGGTTTAAACGTAACTGCGGCCCCGCTTAATTTGGATATGTCGGTGGGTTCCGGCATATGTAGATGCTTGTTGTAACCGAGATGCCTCAGGCAGATACCTTAATGCGACGAAAGGCAGCACTTGTGCTCCATCTAGTTTGACGTATCCCAAGGATGAGATACATACATGAGTGCTCCTCTACTGACACGTTTCGCTTTGCTCACAGCAAAACATTAATCCAACGCAGTCCGCAGGTATGGTGACTAGCGCAAAGTTTGTCTGTATCTTAGTAAGCCGTTAGTTTCGAAGACTGCCGCTACTCTGTTGAACCCATATTCGAACCCTGAAGTCGAAGACGTTCTTGCTCAGCTTGAGAGCCCCTCCGCGCCACGCATCACCCAGTGCCGCTGATGCCCAAGCACAGAAAACGGATGTTGCTATAGAATCCACGGTGTAGGCGAATAATCCATTTTGTCACCCTCAACAACCGACGCTCGTGAGTTCAGGTGAAGTACGGCTTCCTCGTGTTACATACACTTTTACGTATTTGAACTCCGGTATCTACACATTACGAGACGCATTATCAGCGTATCTTGGGCTTTAACGTGTATAGGGCGCCTCAGGATTTGCGGTTATATTTAACGCGCTCTCCTTCCCGCTTCAGGGAATAATAGCAAGCGTGTTTTTTAGGAAAGTCAAATGCATGGATCAGGGGCAAGATGCAGACACGGCTTACTTCCCATGGCAGCTTATTGGGTGGATGCCATTC"""
#reference = reference.lower()
#query = "AGTTTGACGTATCCCAAGGATGAGATACAGTACATGAGTGCTCCTCTAC"
#query = query.lower()
reads = []
#reads.append(query)

reads = read_reads('reads.txt')

#read = "CAAGGTTTAAACGTAACTGCGGCCCCGCTTAATTGGATATGTCGCTGGTT"
#reads.append(read)
#print("reads: ", reads)
reference = read_reference('genome.txt')
reference = reference.lower()
data = BWA(reference)

for query in reads:
	difference_threshold = 0
	query = query.lower()
	matches = []
	#matches = data.find_match(query,0)
	while len(matches) == 0:
		matches = data.find_match(query,difference_threshold)
		difference_threshold+=1
	if len(matches)>1:
		print("length of matches is greater than 1: ", matches)
		matches.append(matches[0])
	print(query)
	print("matches: ", matches)
	print("difference threshold: ", difference_threshold-1)
	findMutations(query, matches[0], difference_threshold-1, 0, data)
	print("found mutations")

print("insertionsDict: ", insertionsDict)
print("deletionsDict: ", deletionsDict)
print("substitutionsDict: ", substitutionsDict)

#if __name__ == "__main__":
#	#index the reference and build up all the necessary datastructures (5 of them; BWT, SA, reference alphabet, C and OCC arrays)
#	data = BWA(reference)
#	print ("\n\nReference: \"%s\"" % reference)
#	if show_data_structures:
#		#printing out the datastructues for manual inspection	
#		print ("\nSA		BWT")
#		print ("--		---")
#		for i in range(len(data.SA)):
#			print ("%s		 %s" % (data.SA[i],data.BWT[i]))
#		print ("\nC(a) = the number of characters 'a' in the Reference that are lexographically smaller than 'a'")
#		print (data.C)
#		print ("\nOcc(a,i) is the number of occurances of the character 'a' in BWT[0,i]")
#		print (data.Occ)
#		print ("\nOcc_reverse(a,i) is the number of occurances of the character 'a' in BWT_reverse[0,i]")
#		print( data.Occ_reverse)
#	if indels_allowed: extra = " and insertions and deletions allowed"
#	else: extra = " with no insertions or deletions allowed"
#	print ("Searching for \"%s\" with max difference threshold of %d%s..." % (query,difference_threshold,extra))
#	matches = []
#	while len(matches) == 0:
#		matches = data.find_match(query,difference_threshold)
#		difference_threshold+=1
#	if len(matches)>1:
#		matches = matches[0]

#	if show_data_structures:
#		print ("D array:")
#		print (data.D)
#	print(len(matches))
#	print(matches)
#	print(difference_threshold-1, " number of mismathces")
#	print ("%d match(es) at position(s): %s \n\n" % (len(matches),matches))

#	findMutations(query, matches[0], difference_threshold-1, 0, data)
