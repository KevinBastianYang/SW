#!/Users/jessie/anaconda3/bin/python
#__author__ = "Jucnen Yang" __email__ = "junchen.yang@yale.edu" __copyright__ = "Copyright 2020" __license__ = "GPL"
#__version__ = "1.0.0"
### Implement your Smith-Waterman Algorithm
def runSW(inputFile, scoreFile, openGap=-2, extGap=-1):
	


	########################################################preprocessing:
	####open the input file
	try:
		with open(inputFile,"r") as f:
			lines = f.readlines()
			if len(lines) > 2:
				raise Exception("Mutiple sequences (>= 3) detected: " + str(len(lines)))
			else:
				seq1 = lines[0].strip()
				seq2 = lines[1].strip()
	except IOError:
		print("Input file can't be opened, check the path or the existence of the file.")
		return(-1)
	####open the score file and store in score_schem
	try:
		with open(scoreFile,'r') as fp:
			names = fp.readline().strip().split()	        
			lines = fp.readlines()
			counter = 0
			scores = []
			for line in lines:
				counter +=1
				score = np.array(line.strip().split()[1:])
				if len(score):
					scores.append(score.astype(np.int))
			scores = np.array(scores)
		score_schem = scores
		####name_index used for amino acid index in the score_schem 
		name_index = {nam: idx for idx, nam in enumerate(names)}
	except IOError:
		print("scoreFile can't be opened, check the path or the existence of the file.")
		return(-2)
	####process the openGap and extGap
	try:
		openGap = float(openGap)
	except ValueError:
		print("openGap can't be converted to float, check the input.")
		return(-3)

	try:
		extGap = float(extGap)
	except ValueError:
		print("extGap can't be converted to float, check the input.")
		return(-3)



	########################################################initialization:
	####initialization of the scoring matrix (n+1,m+1)
	score_mat = np.zeros((len(seq2)+1,len(seq1)+1))
	####initialization the similarity matrix (n+1,m+1)
	similarity_mat = np.zeros((len(seq2)+1,len(seq1)+1))
	for i in range(1,len(seq2)+1):
		for j in range(1,len(seq1)+1):
			similarity_mat[i][j] = score_schem[name_index[seq2[i-1]]][name_index[seq1[j-1]]]	
    	####initialization of the trace matrix (n+1,m+1), in which [x,y] represents the trace back position 
    	####note: [-1, -1] represents the current pos in the score_mat is 0
	trace_mat = np.zeros((len(seq2)+1,len(seq1)+1,2))
	for i in range(len(seq2)+1):
		trace_mat[i][0] = np.array([-1,-1])
	for j in range(1,len(seq1)+1):
		trace_mat[0][j] = np.array([-1,-1])




	########################################################alignment:
	####token for best score and its position
	max_value = 0
	max_pos = [0,0]

	####start aligning: from upper left to bottom right 
	####note: the first row and the first column remains 0
	for i in range(1,len(seq2)+1):
		for j in range(1,len(seq1)+1):
			#### find the horizontal best score source, with index max_h_index
			h_values = [score_mat[i][k] + openGap + (j-k-1)*extGap for k in range(j)]
			max_h_index = np.argmax(h_values)
			#### find the vertical best score source, with index max_v_index
			v_values = [score_mat[k][j] + openGap + (i-k-1)*extGap for k in range(i)]
			max_v_index = np.argmax(v_values)
			####get the best score from the optimal-substructure, remember the direction (score_source)
			substructure = [0,score_mat[i-1][j-1]+similarity_mat[i][j],h_values[max_h_index],v_values[max_v_index]]
			score_source = np.argmax(substructure)
			score_mat[i][j] = substructure[score_source]
			####update the best score
			if score_mat[i][j] >= max_value:
				max_value = score_mat[i][j]
				max_pos = [i,j]
			####store the traceback direction
			trace_map = {0:np.array([-1,-1]),1:np.array([i-1,j-1]),2:np.array([i,max_h_index]),3:np.array([max_v_index,j]) }
			trace_mat[i][j] = trace_map[score_source]
	trace_mat = trace_mat.astype(int)




	########################################################traceback:
	current_pos = [max_pos[0],max_pos[1]]
	prev_pos = [trace_mat[current_pos[0]][current_pos[1]][0],trace_mat[current_pos[0]][current_pos[1]][1]]
	####three alignment output strings
	align_str_1 = ""
	align_str_2 = ""
	middle_slash = ""
	####keep tracing back until 0 score
	while prev_pos[0] != -1 or prev_pos[1] != -1:
		####handle horizontal tracing 
		if prev_pos[0] == current_pos[0]:
			align_str_2 += '-'*(current_pos[1]-prev_pos[1])
			for i in range(1,current_pos[1]-prev_pos[1]+1):
				align_str_1 += seq1[current_pos[1]-i]	
			middle_slash += ' '*(current_pos[1]-prev_pos[1])
		####handle vertical tracing
		elif prev_pos[1] == current_pos[1]:
			align_str_1 += '-'*(current_pos[0]-prev_pos[0])
			for i in range(1,current_pos[0]-prev_pos[0]+1):
				align_str_2 += seq2[current_pos[0]-i]
			middle_slash += ' '*(current_pos[0]-prev_pos[0])
		####handle diagonal tracing
		else:
			align_str_2 += seq2[current_pos[0]-1]
			align_str_1 += seq1[current_pos[1]-1]
			if seq2[current_pos[0]-1] == seq1[current_pos[1]-1]:
				middle_slash += '|'
			else:
				middle_slash += ' '
		current_pos = prev_pos
		prev_pos = [trace_mat[current_pos[0]][current_pos[1]][0],trace_mat[current_pos[0]][current_pos[1]][1]]
	####process the 0 score pos
	if current_pos[0] == 0 and current_pos[1] != 0:
		align_str_2 += '-'
		align_str_1 += seq1[current_pos[1]-1]
		middle_slash += ' '
	elif current_pos[1] == 0 and current_pos[0] != 0:
		align_str_1 += '-'
		align_str_2 += seq2[current_pos[0]-1]
		middle_slash += ' '
	elif current_pos[0] != 0 and current_pos[1] != 0:
		align_str_2 += seq2[current_pos[0]-1]
		align_str_1 += seq1[current_pos[1]-1]
		if seq2[current_pos[0]-1] == seq1[current_pos[1]-1]:
				middle_slash += '|'
		else:
			middle_slash += ' '


	#with open("SW_align_results.txt","ab") as fw:
		#f.write("----------------------------------------------------\n")

	########################################################output:
	try:
		with open("SW_align_results.txt","a") as fw:
			fw.write("----------------------------------------------------\n")
			fw.write("|		Sequences                         |\n")
			fw.write("----------------------------------------------------\n")
			fw.write('\n')
			fw.write("seq1\n")
			fw.write(seq1+'\n')
			fw.write("seq2\n")
			fw.write(seq2+'\n')
			fw.write("----------------------------------------------------\n")
			fw.write("|		Score Matrix                      |\n")
			fw.write("----------------------------------------------------\n")	
			fw.write('\t\t')
			for i in list(seq1):
				fw.write(i+'\t')
			fw.write('\n')
			fw.write('\t')
			for i in range(len(seq2)):
				for j in range(len(seq1)+1):
					fw.write(str(score_mat[i][j]))
					fw.write('\t')
				fw.write('\n')
				fw.write(list(seq2)[i]+'\t')
			for j in range(len(seq1)+1):
				fw.write(str(score_mat[len(seq2)][j]))
				fw.write('\t')
			fw.write('\n')


			fw.write("----------------------------------------------------\n")
			fw.write("|        Best Local Alignment                      |\n")
			fw.write("----------------------------------------------------\n")

			fw.write("Alignment Score:")
			fw.write(str(max_value)+'\n')
			fw.write("Alignment Results:\n")
			fw.write(align_str_1[::-1]+'\n')
			fw.write(middle_slash[::-1]+'\n')
			fw.write(align_str_2[::-1]+'\n')		

	except IOError:
		print("output file can't be opened.")
		return(-4)

	print("Aligment finished! Alignment results are in file: SW_align_results.txt\n")
	return(0)
