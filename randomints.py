from random import random

f = open('inputfile', 'w')
for x in xrange(100):
    y = int(random()*(10**12)) + 1
    s = str(y)
    f.write(s + "\n")