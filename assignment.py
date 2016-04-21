__author__ = 'nishaswarup'

import structures
from random import random
from random import sample
from time import time

'''
The order of this file is
    - KK Algorithm
    - Functions for Representation 1
    - Functions for Representation 2
    - Functions for the 3 algorithms
    - Functions for Testing
'''

def KK(nums):
    '''
    This function runs the KK algorithm on the list nums
    It builds a max-heap from the numbers, pops the largest 2
    numbers, then pushes their difference. It does this until
    the heap has size 1, and then it returns that value.
    This is O(n*log(n)) time. Building the heap is a lot of pushes,
    which is log(n) for each heap, so n*log(n) total. Then, popping
    is log(n, and so is pushing. We do each of those 3N times, so we
    get another O(n*log(n)). So the overall algorithm is O(n*log(n))
    :param nums: a list of numbers
    :return: A single value, which is the result of the KK algorithm
    '''

    # first we construct the heap
    heap = []
    for x in nums:
        structures.heap_push(heap, x)
    # now we've constructed the heap, now we loop through
    # and replace the largest two elements with their difference
    while len(heap) > 1:
        biggest = structures.heap_pop(heap)
        biggest2 = structures.heap_pop(heap)
        structures.heap_push(heap, biggest - biggest2)
    return heap[0]

'''
    ****** Representation 1 *******
    *  The below functions are all for the first representation, where
    *  we generate solutions that are alternating -1 and 1
    *******************************
'''

def calculate_residue_1(nums, solution):
    '''
    This function calculates the residue for a given solution
    :param nums: A list of numbers that we're trying to find the parition of
    :param solution: A list of -1 and 1
    :return: A number which represents the residue
    '''
    sum = 0
    for index, value in enumerate(nums):
        sum += value*solution[index]
    return abs(sum)

def generate_random_soln_1(nums):
    '''
    This function generates a random solution for the first representation
    So, it just returns a list of -1 and 1 that is the same length as nums
    :param nums: A list of numbers
    :return: A list of -1 and 1
    '''
    out = []
    for x in xrange(len(nums)):
        if random() < .5:
            out.append(-1)
        else:
            out.append(1)
    return out

def random_move_1(soln):
    '''
    This function makes a random move on the given solution. It retuns
    an array of swaps to make in this random move. For example, it could
    return [[5, 1, -1], [4, -1, 1]]
    Which would indicate we should change position 5 from 1 to -1 and position 4 from -1 to 1.
    The reason we use this strange notation is to make it integrate  with random_move_2
    :param soln: a list of -1 or 1
    :return: An array, look above
    '''
    # Generating two random integers within a range 
    try:
        [i,j] = sample(range(0, len(soln)), 2)
    except ValueError:
        print('Sample size exceeded population size.')


    # This was the code you had, but I don't know if it actually generates two random integers
    # within a range, I found this function and it also checks that we have enough space in the
    # range   TODO: eliminate once you agree

    # i = int(random()*len(soln))
    # j = int(random()*len(soln) - 1)
    # if j >= i: j += 1


    # now i and j are two random numbers s.t. i != j
    out = [[i, soln[i], soln[i]*-1]]
    if random() < .5:
        out.append([j, soln[j], soln[j]*-1])
    return out
'''
    ****** Representation 2 *******
    *  The below functions are all for the second representation, where
    *  we generate solutions that group numbers together
    *******************************
'''

def calculate_residue_2(nums, soln):
    '''
    This calculates the residue for a solution in representation 2.
    :param nums: A list of numbers
    :param soln: A solution in representation 2, sequence representing a partitioning
    :return: A number that represents the residue
    '''

    for i in range(len(soln)):
        if nums[i] != 0:
            for j in range(i+1,len(soln)):
                if soln[i] == soln[j]:
                    nums[i] += nums[j]
                    nums[j] = 0
    return KK(nums)







def generate_random_soln_2(nums):
    '''
    We generate a random solution for the second representation.
    This is just a list of numbers from 0 to len(nums)
    :param nums: A list of numbers
    :return: A list of numbers of the same size as nums
    '''
    out = []
    for x in xrange(nums):
        out.append(int(random()*len(nums)))
    return out

'''
    ****** 3 algorithms *******
    *  The below functions are the 3 algorithms we needed to implement
    ***************************
'''
def repeated_random(nums, max_iter, generate_soln, find_residue):
    '''
    This function runs the first algorithm given, where we
    find lots of random solutions and take the best one
    :param nums: The list of numbers we're finding the number partition for
    :param max_iter: The number we're iterating to
    :param generate_soln: A function for generating random solutions (this varies for represenations 1 and 2)
    :param find_residue: A function for calculating the residue of a solution (again, varies per representation)
    :return: Soln, residue. Solution is a list of numbers, residue is the best residue
    '''
    # generate a random solution
    best_soln = generate_soln(nums)
    best_residue = find_residue(nums, best_soln)
    # we loop through lots of times, and return the
    # best random solution
    for x in xrange(max_iter):
        soln = generate_soln(nums)
        residue = find_residue(nums, soln)
        if residue < best_residue:
            best_residue = residue
            best_soln = soln
    return best_soln, best_residue

def hill_climb(nums, max_iter, generate_soln, find_residue, find_neighbor):
    '''
    This function runs the second algorithm, where we find a random solution,
    and then make local improvements.
    :param nums: The list of numbers we're trying to find the partition for
    :param max_iter: The number we're iterating to
    :param generate_soln: A function used to generate a solution
    :param find_residue: A function used to calculate a residue
    :param find_neighbor: A function that finds a random neighbor. Should be either random_move_1 or random_move_2
    :return: Soln, residue. Solution is a list of numbers, residue is the best residue
    '''
    best_soln = generate_soln(nums)
    best_residue = find_residue(nums, best_soln)
    for x in xrange(max_iter):
        # find the moves we should make to get to a random neighbor
        moves = find_neighbor(best_soln)
        # make these moves
        for move in moves:
            best_soln[move[0]] = move[2]
        residue = find_residue(nums, best_soln)
        # if this was a good move, we keep it
        if residue < best_residue:
            best_residue = residue
        # otherwise we revert
        else:
            for move in moves:
                best_soln[move[0]] = move[1]
    return best_soln, best_residue

'''
    ****** Testing *******
    *  The below functions are all for testing
    **********************
'''
def generate_random_instance():
    '''
    This returns a list of 100 numbers, with each number being chosen from [1, 10**12]
    This specification is given in the assignment
    :return: A list of numbers
    '''
    out = []
    for x in xrange(100):
        out.append( int(random()*(10**12)) + 1 )
    return out

def test_instance(max_iter):
    '''
    This function tests one instance with the KK algorithm, a repeated random algorithm,
    a hill climbing algorithm, and a simulated annealing algorithm for 25,000 iterations.
    :param max_iter: how many iterations we should run
    :return: Notihg, for now
    '''
    # generate instance
    instance = generate_random_instance()
    # KK algorithm
    print "KK Algorithm"
    start_time = time()
    KK_ans = KK(instance)
    KK_time = time() - start_time
    print "\t\tResult:\t" + str(KK_ans)
    print "\t\tTime:\t" + str(KK_time) + "s"
    # Run Representation 1 Algorithm
    print "Representation 1:"
    # repeated random
    print "\tRepeated Random"
    start_time = time()
    temp,R1_RR_ans = repeated_random(instance, max_iter, generate_random_soln_1, calculate_residue_1)
    R1_RR_time = time() - start_time
    print "\t\tResult:\t" + str(R1_RR_ans)
    print "\t\tTime:\t" + str(R1_RR_time)
    # hill climbing
    print "\tHill Climbing"
    start_time = time()
    temp, R1_HC_ans = hill_climb(instance, max_iter, generate_random_soln_1, calculate_residue_1, random_move_1)
    R1_HC_time = time() - start_time
    print "\t\tResult:\t" + str(R1_HC_ans)
    print "\t\tTime:\t" + str(R1_HC_time) + "s"

def test_calculate_residue_2 ():
    '''
    This function tests a small number of cases for which the solution has been 
    calculated by hand
    :param: nothing for now
    :return: boolean, true if it passes the tests
    '''
    nums = [10,8,7,6,5]
    soln = [1,2,2,4,5]
    return 4 == calculate_residue_2(nums,soln)
    
test_instance(25000)
test_calculate_residue_2()
