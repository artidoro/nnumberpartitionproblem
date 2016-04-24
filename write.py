file = open('output.txt','w')
for i in range(99):
	file.write(str(i)+'\n')
file.write('99')
file.close()