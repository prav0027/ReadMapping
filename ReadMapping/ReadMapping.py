import pandas as pd

class BWA:
	#this initializer function creates all the datastructures necessary using the reference string
	def __init__(self, reference):
		#declare datastructures
		rotation_list = list()
		rotation_list_reverse = list()
		suffix_array = list()
		bwt = list()
		C = dict()
		Occ = dict()
		Occ_reverse = dict()
		alphabet = set()
		reverse_reference = reference[::-1]#reverse reference
		
		alphabet.add('c')
		alphabet.add('g')
		alphabet.add('t')
		alphabet.add('a')
		
		for char in alphabet:
			C[char] = 0
			Occ[char] = list()# from website: in Occ, each character has an associated list of integer values
			Occ_reverse[char] = list()
	
		#add ending character to reference
		reference = reference+"$"
		reverse_reference = reverse_reference+"$"

		#create suffix combinations of the reference and reverse reference
		for i in range(len(reference)):
			new_rotation = "%s%s" % (reference[i:],reference[0:i])
			suffixObj = Suffix(new_rotation,i)
			rotation_list.append(suffixObj)
			
			new_rotation_reverse = "%s%s" % (reverse_reference[i:],reverse_reference[0:i])
			suffixObj_rev = Suffix(new_rotation_reverse,i)
			rotation_list_reverse.append(suffixObj_rev)
		
			#number of characters that are smaller than each char
			if reference[i]!='$':
				for char in alphabet:
					if reference[i] < char:
						C[char] = C[char] + 1	
		
		#sort the suffixes
		rotation_list.sort(key=textKey)
		rotation_list_reverse.sort(key=textKey)
	
		#record results into the suffix array and the BWT array and calculate Occ
		for i in rotation_list:
			suffix_array.append(i.pos)
			bwt.append(i.text[-1:])
		
			#create the Occ
			for char in alphabet:
				if len(Occ[char]) == 0:
					prev = 0
				else:
					prev = Occ[char][-1]
				if i.text[-1:] == char:
					Occ[char].append(prev+1)
				else:
					Occ[char].append(prev)
					
		#record results into the suffix array and the BWT array and calculate Occ
		for i in rotation_list_reverse:
			#create the Occ
			for char in alphabet:
				if len(Occ_reverse[char]) == 0:
					prev = 0
				else:
					prev = Occ_reverse[char][-1]
				if i.text[-1:] == char:
					Occ_reverse[char].append(prev+1)
				else:
					Occ_reverse[char].append(prev)					
					
		self.SA = suffix_array
		self.BWT = bwt
		self.C = C
		self.Occ = Occ
		self.Occ_reverse = Occ_reverse
		self.n = len(reference)
		self.D = list()#empty list for later use
		self.alphabet = alphabet

	#get position of the query in the reference
	def find_match(self,query,num_differences):
		if num_differences == 0:
			return self.exact_match(query)
		else:
			return self.inexact_match(query,num_differences)


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

	#threshold is number of allowed mutations
	def inexact_match(self,query,threshold):
		self.calculate_d(query)
		SA_indeces = self.inexact_recursion(query, len(query)-1, threshold, 0, self.n-1)
		return [self.SA[x] for x in SA_indeces]

	def inexact_recursion(self,query,i,threshold,k,l):
		tempset = set()
		resultZ = dict()
		resultDict = dict()
		if (threshold < self.get_D(i)):
			return set()
		if i < 0:
			for m in range(k,l+1):
				tempset.add(m)
				resultZ[m] = threshold
				resultDict[m] = resultZ
			return tempset
			
		result = set()
		#insertion
		resultAdd = self.inexact_recursion(query,i-1,threshold-insertion_penalty,k,l)
		result = result.union(resultAdd)
		for char in self.alphabet:
			newK = self.C[char] + self.OCC(char,k-1) + 1 
			newL = self.C[char] + self.OCC(char,l)
			if newK <= newL:#if the substring was found
				resultAdd = self.inexact_recursion(query,i,threshold-deletion_penalty,newK,newL)
				result = result.union(resultAdd)
				if char == query[i]:#if the char was correctly aligned, then continue without decrementing z (differences)
					result = result.union(self.inexact_recursion(query,i-1,threshold,newK,newL))
				else:#continue but decrement z, to indicate that this was a difference/unalignment
					result = result.union(self.inexact_recursion(query,i-1,threshold-mismatch_penalty,newK,newL))
		return result

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
insertionsDictNoErrors = dict()
deletionsDict = dict()
deletionsDictNoErrors = dict()
substitutionsDict = dict()
substitutionsDictNoErrors = dict()

matchIndexes = list()

def findMutations(query, index, numberOfMismatches, mismatchesCount, data):
	#mismatchesCount=0
	for i in range(len(query)):
		substitutionNeeded = True
		mismatchesCountTemp = mismatchesCount
		insertions_i_list = []
		insertions_char_list = []
		deletions_i_list = []
		if reference[index]==query[i]:
			index+=1
			continue
		#check for insertions
		while numberOfMismatches>0:
			iTemp = i+1
			if index in data.find_match(query[iTemp:], numberOfMismatches-1):
				print("insertion @ index: ", index)
				substitutionNeeded = False
				insertions_char_list.append(query[i])
				insertions_i_list.append(index-1)
				index+=-1
				numberOfMismatches+=-1
				mismatchesCountTemp+=1
				mismatchesCount = mismatchesCountTemp
			else:
				break
		for i in range(len(insertions_i_list)):
			if insertions_i_list[i] not in list(insertionsDict.keys()):
				insertionsDict[insertions_i_list[i]] = list(insertions_char_list[i])
			else:
				tempList = insertionsDict[insertions_i_list[i]]
				tempList.append(insertions_char_list[i])
				insertionsDict[insertions_i_list[i]] = tempList
		while numberOfMismatches>0:
			if index+1 in data.find_match(query[i:], numberOfMismatches-1):
				substitutionNeeded = False
				deletions_i_list.append(index)
				index+=1
				numberOfMismatches+=-1
				mismatchesCountTemp+=1
				mismatchesCount = mismatchesCountTemp
			else:
				break
		for val in deletions_i_list:
			if val not in list(deletionsDict.keys()):
				deletionsDict[val] = list(reference[val])
			else:
				tempList = deletionsDict[val]
				tempList.append(reference[val])
				deletionsDict[val] = tempList

		if substitutionNeeded==True:
			nucleotide_tuple = (reference[index], query[i])
			numberOfMismatches+=-1
			if index not in list(substitutionsDict.keys()):
				substitutionTempList = []
				substitutionTempList.append(nucleotide_tuple)
				substitutionsDict[index] = substitutionTempList
			else:
				tempList = substitutionsDict[index]
				tempList.append(nucleotide_tuple)
				substitutionsDict[index] = tempList
		index+=1
	return

def takeOutErrorReads():
	deletionsDictTemp = dict()
	insertionsDictTemp = dict()
	substitutionsDictTemp = dict()
	for i in list(deletionsDict.keys()):
		minVal = max(0,i-50)
		deletionsCount = 0
		for k in matchIndexes:
			if k >= minVal and k <= i:
				deletionsCount+=1
		deletionsPercent = len(deletionsDict[i])/deletionsCount
		if deletionsPercent>0.5:
			tempList = deletionsDict[i]
			deletionsDictTemp[i] = tempList[0]
	#print("deletions temp: ", deletionsDictTemp)
	myKeys = list(deletionsDictTemp.keys())
	myKeys.sort()
	#print(myKeys)
	deletionsNoErrors = {i: deletionsDictTemp[i] for i in myKeys}
	#print("deletionsDictNoErrors: ", deletionsNoErrors)
	for i in list(insertionsDict.keys()):
		minVal = max(0, i-50)
		insertionsCount = 0
		for k in matchIndexes:
			if k >= minVal and k<=i:
				insertionsCount+=1
		insertionsPercent = len(insertionsDict[i])/insertionsCount
		if insertionsPercent>0.5:
			#check all of the elements in the list are the same!
			tempList = insertionsDict[i]
			insertionsDictTemp[i] = tempList[0]
	#print("insertions temp: ", insertionsDictTemp)
	myKeys = list(insertionsDictTemp.keys())
	myKeys.sort()
	#print(myKeys)
	insertionsNoErrors = {i:insertionsDictTemp[i] for i in myKeys}
	#print("insertionsDictNoErrors: ", insertionsDictNoErrors)
	for i in list(substitutionsDict.keys()):
		minVal = max(0,i-50)
		substitutionsCount=0
		for k in matchIndexes:
			if k >= minVal and k<=i:
				substitutionsCount+=1
		substitutionsPercent = len(substitutionsDict[i])/substitutionsCount
		if substitutionsPercent>0.5:
			#check all of the elements in the list are the same!
			tempList = substitutionsDict[i]
			substitutionsDictTemp[i] = tempList[0]	
	#print("substitutions temp: ", substitutionsDictTemp)
	myKeys = list(substitutionsDictTemp.keys())
	myKeys.sort()
	#print(myKeys)
	substitutionsNoErrors = {i:substitutionsDictTemp[i] for i in myKeys}
	#print("substitutionsDictNoErrors: ", substitutionsDictNoErrors)
	return deletionsNoErrors, insertionsNoErrors, substitutionsNoErrors

def save_results(fileName = 'results_2.csv'):
	f = open(fileName, 'w')
	for i in list(substitutionsDictNoErrors.keys()):
		tempTuple = substitutionsDictNoErrors[i]
		f.write('>S'+str(i)+' '+tempTuple[0].upper()+' '+tempTuple[1].upper()+'\n')
	for i in list(insertionsDictNoErrors.keys()):
		f.write('>I'+str(i)+' '+insertionsDictNoErrors[i].upper()+'\n')
	for i in list(deletionsDictNoErrors.keys()):
		f.write('>D'+str(i)+' '+deletionsDictNoErrors[i].upper()+'\n')
	f.close()

	#read_file = 

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

insertion_penalty = 1
deletion_penalty = 1
mismatch_penalty = 1

reads = []

reads = read_reads('reads.txt')


reference = read_reference('genome.txt')
reference = reference.lower()
data = BWA(reference)

count = 0
for query in reads:
	difference_threshold = 0
	query = query.lower()
	matches = []
	while len(matches) == 0:
		matches = data.find_match(query,difference_threshold)
		difference_threshold+=1
	if len(matches)>1:
		#print("length of matches is greater than 1: ", matches)
		matches.append(matches[0])

	matchIndexes.append(matches[0])

	#print(query)
	#print("matches: ", matches)
	#print("difference threshold: ", difference_threshold-1)
	findMutations(query, matches[0], difference_threshold-1, 0, data)
	print("count: ", count)
	count+=1
	#print("found mutations")

deletionsDictNoErrors, insertionsDictNoErrors, substitutionsDictNoErrors = takeOutErrorReads()
save_results()

print("insertionsDict: ", insertionsDictNoErrors)
print("deletionsDict: ", deletionsDictNoErrors)
print("substitutionsDict: ", substitutionsDictNoErrors)
