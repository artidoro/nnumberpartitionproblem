import structures
from random import random
from random import sample
from time import time


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
    # try:
    #     [i,j] = sample(range(1, len(soln)+1), 2)
    # except ValueError:
    #     print('Sample size exceeded population size.')


    # This was the code you had, but I don't know if it actually generates two random integers
    # within a range, I found this function and it also checks that we have enough space in the
    # range   TODO: eliminate once you agree

    # i = int(random()*len(soln))
    # j = int(random()*len(soln) - 1)
    # if j >= i: j += 1

    i,j = 0,0
    [i,j] = sample(range(0, len(soln)), 2)
    # now i and j are two random numbers s.t. i != j
    out = [[i, soln[i], soln[i]*-1]]
    if random() < .5:
        out.append([j, soln[j], soln[j]*-1])
    return out


print random_move_1([1,1,1,1])