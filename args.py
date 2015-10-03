i = 0
mypi = 0.0
while round(mypi*4,2) != 3.14:
    mypi += ( (-1)**i / (2.0*i + 1) )
    i += 1
mypi *= 4
print "using a sum of %d terms my estimate of pi is rounded to 3.14" % i