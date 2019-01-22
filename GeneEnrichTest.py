import scipy.stats as st


# total number of genes in the list
ar1 = 2000

# total number of genes in the list with condition A
ar2 = 500
 

# total number of genes in the list with condition B
ar3 = 100


# number of genes with both condition A and B
ar4 = 70

#~ print 'total number in population: ' + sys.argv[1]
#~ print 'total number with condition in population: ' + sys.argv[2]
#~ print 'number in subset: ' + sys.argv[3]
#~ print 'number with condition in subset: ' + sys.argv[4]
p_val_minus = st.hypergeom.cdf(int(ar4),int(ar1),int(ar2),int(ar3))
p_val_plus  = st.hypergeom.sf(int(ar4) - 1,int(ar1),int(ar2),int(ar3))



#~ The result will be

#~ a p-value where by random chance number of genes with both condition A and B will be <= to your number with condition A and B
#~ a p-value where by random chance number of genes with both condition A and B will be >= to your number with condition A and B
#~ The second p-value is probably what you want.



def pValue(tot, num_a, num_b, num_ab):
    return st.hypergeom.sf(int(num_ab) - 1,int(tot),int(num_a),int(num_b))



print 'p-value <= ' + str(p_val_minus)
print 
print 'p-value >= ' + str(p_val_plus)

print 'p-value >= ' + str(pValue(ar1,ar2,ar3,ar4))
print 'p-value >= ' + str(pValue(ar1,ar3,ar2,ar4))

