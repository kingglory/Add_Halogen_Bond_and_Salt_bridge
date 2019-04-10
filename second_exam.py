import string
def s():
  str = raw_input("string here:")
  xiaokuohao = "()"
  zhongkuohao = "[]"
  dakuohoa = "{}"
  x0=[]
  x1=[]
  z0=[]
  z1=[]
  d0=[]
  d1=[]
  for i,s in enumerate(str):
    if s == xiaokuohao[0]:
      x0.append((i,s))
    if s == xiaokuohao[1]:
      x1.append((i,s))
    if s == dakuohoa[0]:
      d0.append((i, s))
    if s == dakuohoa[1]:
      d1.append((i, s))
    if s == zhongkuohao[0]:
      z0.append((i, s))
    if s == zhongkuohao[1]:
      z1.append((i, s))
  for i,s in enumerate(str):
    for j,k in enumerate(x0):
      if s == k[1]:
        n1 = k[0]
        n2 = x1[j][0]
        string = str[n1+1:n2]
        string1 = str[0:n1]
        l1 = len(string1)
        for h in range(l1):
          if string1[h:l1].isdigit():
            number = int(string1[h:l1])
            strs = str.replace(str[n1:n2 + 1], (number) * (string))


  for i,s in enumerate(str):
    for j,k in enumerate(z0):
      if s == k[1]:
        n1 = k[0]
        n2 = z1[j][0]
        string = str[n1+1:n2]
        string1 = str[0:n1]
        l1 = len(string1)
        for h in range(l1):
          if string1[h:l1].isdigit():
            number = int(string1[h:l1])
            strs = str.replace(str[n1:n2 + 1], (number) * (string))

  for i,s in enumerate(str):
    for j,k in enumerate(d0):
      if s == k[1]:
        n1 = k[0]
        n2 = d1[j][0]
        string = str[n1+1:n2]
        string1 = str[0:n1]
        l1 = len(string1)
        for h in range(l1):
          if string1[h:l1].isdigit():
            number = int(string1[h:l1])
            strs = str.replace(str[n1:n2 + 1], (number) * (string))

  return strs
strs = s()
l = list(strs)
new_string = l.reverse()
result = ''.join(l)

print result