infile = open('input.txt')
soln = []
for line in infile:
	soln.append(line)
soln = [int(i) for i in soln]
print soln
infile.close()	
