# Import the packages
import numpy as np
from scipy import stats


## Define 2 random distributions


#Gaussian distributed data with mean = 2 and var = 1
#a = np.random.randn(N) + 2
a = np.array([41.82338062, 42.45045037, 55.18443692])


#Gaussian distributed data with with mean = 0 and var = 1
#b = np.random.randn(N)
b = np.array([10178.90062826, 12275.42266438,  9697.02723771])
#Sample Size
N = len(a)



## Calculate the Standard Deviation
#Calculate the variance to get the standard deviation

#For unbiased max likelihood estimate we have to divide the var by N-1, and therefore the parameter ddof = 1
var_a = a.var(ddof=1)
var_b = b.var(ddof=1)

#std deviation
s = np.sqrt((var_a + var_b)/2)



## Calculate the t-statistics
t = (a.mean() - b.mean())/(s*np.sqrt(2/N))



## Compare with the critical t-value
#Degrees of freedom
df = 2*N - 2

#p-value after comparison with the t 
p =  stats.t.cdf(-abs(t),df=df)

p11 = 2*(1-p)
p12 = 2*p

print("t = " + str(t))
print("p0 = " + str(p))
print("p11 = " + str(p11))
print("p12 = " + str(p12))
#Note that we multiply the p value by 2 because its a twp tail t-test
### You can see that after comparing the t statistic with the critical t value (computed internally) we get a good p value of 0.0005 and thus we reject the null hypothesis and thus it proves that the mean of the two distributions are different and statistically significant.


## Cross Checking with the internal scipy function
t2, p2 = stats.ttest_ind(a,b)
print("t = " + str(t2))
print("p = " + str(p2))
