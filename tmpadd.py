import sys

a = b = c = 0
with open(sys.argv[1]) as f:
   for line in f:
      items = line.split()
      a += int(items[2])
      b += int(items[3])
      c += int(items[4])
print a, b, c
