import structures
from random import random
from time import time




def residue (nums, soln):
	for i in range(len(soln)):
		if nums[i] != 0:
			for j in range(i+1,len(soln)):
				if soln[i] == soln[j]:
					nums[i] += nums[j]
					nums[j] = 0
	return nums


nums = [1,1,1,1,1,1,1,1,1,1]
soln1 = [1,1,3,4,5,6,7,8,9,10]
soln2 = [1,1,1,4,1,1,7,8,1,10]
soln3 = [1,1,3,4,10,6,7,8,9,10]
soln4 = [1,1,3,4,5,1,7,8,9,1]


print residue([1,1,1,1,1,1,1,1,1,1], range(1,11))
print residue([1,1,1,1,1,1,1,1,1,1], soln1)
print residue([1,1,1,1,1,1,1,1,1,1], soln2)
print residue([1,1,1,1,1,1,1,1,1,1], soln3)
print residue([1,1,1,1,1,1,1,1,1,1], soln4)