import pandas as pd

class BWT:
	def __init__(self, reference):
		suffixRotation = list()
		reverseSuffixRotation = list()
		suffix_array = list()
		bwt = list()
		C = dict()
		Occ = dict()
		Occ_reverse = dict()
		alphabet = set()
		referenceReverse = reference[::-1]
		
		alphabet.add('c')
		alphabet.add('g')
		alphabet.add('t')
		alphabet.add('a')
		
		for letter in alphabet:
			C[letter] = 0
			Occ[letter] = list()# from website: in Occ, each character has an associated list of integer values
			Occ_reverse[letter] = list()
	
		#add ending character to reference
		reference = reference+"$"
		referenceReverse = referenceReverse+"$"

		#create suffix combinations of the reference and reverse reference
		for i in range(len(reference)):
			new_rotation = "%s%s" % (reference[i:],reference[0:i])
			suffixObj = Suffix(new_rotation,i)
			suffixRotation.append(suffixObj)
			
			new_rotation_reverse = "%s%s" % (referenceReverse[i:],referenceReverse[0:i])
			suffixObj_rev = Suffix(new_rotation_reverse,i)
			reverseSuffixRotation.append(suffixObj_rev)
		
			#number of characters that are smaller than each char
			if reference[i]!='$':
				for char in alphabet:
					if reference[i] < char:
						C[char] = C[char] + 1	
		
		suffixRotation.sort(key=textKey)
		reverseSuffixRotation.sort(key=textKey)
	
		for i in suffixRotation:
			suffix_array.append(i.pos)
			bwt.append(i.text[-1:])
		
			#create the Occ
			for letter in alphabet:
				if len(Occ[letter]) == 0:
					prev = 0
				else:
					prev = Occ[letter][-1]
				if i.text[-1:] == letter:
					Occ[letter].append(prev+1)
				else:
					Occ[letter].append(prev)
					
		#record results into the suffix array and the BWT array and calculate Occ
		for i in reverseSuffixRotation:
			#create the Occ
			for letter in alphabet:
				if len(Occ_reverse[letter]) == 0:
					prev = 0
				else:
					prev = Occ_reverse[letter][-1]
				if i.text[-1:] == letter:
					Occ_reverse[letter].append(prev+1)
				else:
					Occ_reverse[letter].append(prev)					
					
		self.SA = suffix_array
		self.BWT = bwt
		self.C = C
		self.Occ = Occ
		self.Occ_reverse = Occ_reverse
		self.n = len(reference)
		self.newArray = list()
		self.alphabet = alphabet

	def getMatchIndex(self,read,num_differences):
		if num_differences == 0:
			return self.matchNoDifferences(read)
		else:
			return self.matchWithDifferences(read,num_differences)


	def matchNoDifferences(self, read):
		read = read.lower()
		i = 0
		j = self.n - 1
		
		for x in range(len(read)):
			newChar = read[-x-1]
			newI = self.C[newChar] + self.OCC(newChar,i-1) + 1
			newJ = self.C[newChar] + self.OCC(newChar,j)
			i = newI
			j = newJ
		matches = self.SA[i:j+1]
		return matches

	#threshold is number of allowed mutations
	def matchWithDifferences(self,read,threshold):
		self.calculate_distance(read)
		SA_indeces = self.recursionMatchingWithDifferences(read, len(read)-1, threshold, 0, self.n-1)
		return [self.SA[x] for x in SA_indeces]

	def recursionMatchingWithDifferences(self,read,i,threshold,m,n):
		tempset = set()
		resultZ = dict()
		resultDict = dict()
		if (threshold < self.get_D(i)):
			return set()
		if i < 0:
			for j in range(m,n+1):
				tempset.add(j)
				resultZ[j] = threshold
				resultDict[j] = resultZ
			return tempset
			
		result = set()
		#insertion
		resultAdd = self.recursionMatchingWithDifferences(read,i-1,threshold-insertion_penalty,m,n)
		result = result.union(resultAdd)
		for char in self.alphabet:
			newM = self.C[char] + self.OCC(char,m-1) + 1 
			newN = self.C[char] + self.OCC(char,n)
			if newM <= newN:#if the substring was found
				resultAdd = self.recursionMatchingWithDifferences(read,i,threshold-deletion_penalty,newM,newN)
				result = result.union(resultAdd)
				if char == read[i]:#char aligned correctly
					resultAdd = self.recursionMatchingWithDifferences(read,i-1,threshold,newM,newN)
					result = result.union(resultAdd)
				else:#unaligned
					resultAdd = self.recursionMatchingWithDifferences(read,i-1,threshold-mismatch_penalty,newM,newN)
					result = result.union(resultAdd)
		return result

	def calculate_distance(self,read):
		m = 0
		l = self.n-1
		distance = 0
		self.newArray = list()
		for i in range(len(read)):
			m = self.C[read[i]] + self.OCC(read[i],m-1,reverse=True) + 1
			l = self.C[read[i]] + self.OCC(read[i],l,reverse=True)
			if m > l:
				m = 0
				l = self.n - 1
				distance = distance + 1
			self.newArray.append(distance)

	def OCC(self,char,index,reverse=False):
		if index < 0:
			return 0
		else:
			if reverse:
				return self.Occ_reverse[char][index]
			else:
				return self.Occ[char][index]
	
	def get_D(self,index):
		if index < 0:
			return 0
		else:
			return self.newArray[index]

class Suffix:
	def __init__(self, text, position):
		self.text = text
		self.pos = position

def textKey( a ): return a.text

insertionsDict = dict()
insertionsDictNoErrors = dict()
deletionsDict = dict()
deletionsDictNoErrors = dict()
substitutionsDict = dict()
substitutionsDictNoErrors = dict()

matchIndexes = list()

def findMutations(query, index, numberOfMismatches, mismatchesCount, data):
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
			if index in data.getMatchIndex(query[iTemp:], numberOfMismatches-1):
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
			if index+1 in data.getMatchIndex(query[i:], numberOfMismatches-1):
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
	myKeys = list(deletionsDictTemp.keys())
	myKeys.sort()
	deletionsNoErrors = {i: deletionsDictTemp[i] for i in myKeys}
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
	myKeys = list(insertionsDictTemp.keys())
	myKeys.sort()
	insertionsNoErrors = {i:insertionsDictTemp[i] for i in myKeys}
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
	myKeys = list(substitutionsDictTemp.keys())
	myKeys.sort()
	substitutionsNoErrors = {i:substitutionsDictTemp[i] for i in myKeys}
	return deletionsNoErrors, insertionsNoErrors, substitutionsNoErrors

def save_results(fileName = 'results.txt'):
	f = open(fileName, 'w')
	for i in list(substitutionsDictNoErrors.keys()):
		tempTuple = substitutionsDictNoErrors[i]
		f.write('>S'+str(i)+' '+tempTuple[0].upper()+' '+tempTuple[1].upper()+'\n')
	for i in list(insertionsDictNoErrors.keys()):
		f.write('>I'+str(i)+' '+insertionsDictNoErrors[i].upper()+'\n')
	for i in list(deletionsDictNoErrors.keys()):
		f.write('>D'+str(i)+' '+deletionsDictNoErrors[i].upper()+'\n')
	f.close()

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
data = BWT(reference)

count = 0
for query in reads:
	difference_threshold = 0
	query = query.lower()
	matches = []
	while len(matches) == 0:
		matches = data.getMatchIndex(query,difference_threshold)
		difference_threshold+=1
	if len(matches)>1:
		matches.append(matches[0])

	matchIndexes.append(matches[0])

	findMutations(query, matches[0], difference_threshold-1, 0, data)
	print("count: ", count)
	count+=1

deletionsDictNoErrors, insertionsDictNoErrors, substitutionsDictNoErrors = takeOutErrorReads()
save_results()

print("insertionsDict: ", insertionsDictNoErrors)
print("deletionsDict: ", deletionsDictNoErrors)
print("substitutionsDict: ", substitutionsDictNoErrors)
