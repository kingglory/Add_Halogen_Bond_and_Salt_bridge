#wensong
#2018.10.6
# compute the distance between two atom
import math

def distance(atom1,atom2):
    distance = math.sqrt((atom1[0]-atom2[0])**2 + (atom1[1]-atom2[1])**2 + (atom1[2]-atom2[2])**2)
    return distance


# test

atom1 = [1,1,1]
atom2 = [2,1,1]
d = distance(atom1,atom2)
print d