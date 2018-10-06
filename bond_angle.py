#wensong
#2018.10.6
#to calculate the angle between the two bond

import math

def angle(atom1,atom2,atom3):
    #Calculate the three sides length
    x = math.sqrt((atom2[0]-atom1[0])**2+(atom2[1]-atom1[1])**2+(atom2[2]-atom1[2])**2)
    y = math.sqrt((atom3[0]-atom1[0])**2+(atom3[1]-atom1[1])**2+(atom3[2]-atom1[2])**2)
    z = math.sqrt((atom3[0]-atom2[0])**2+(atom3[1]-atom2[1])**2+(atom3[2]-atom2[2])**2)

    angle = math.degrees(math.acos((x**2+y**2-z**2)/(2*x*y)))
    return angle




#test
atom1 = [0,0,0]
atom2 = [1,1,0]
atom3 = [-1,1,0]
print angle(atom1,atom2,atom3)