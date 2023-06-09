class Gaussian:
    
    def __init__(self, a=0, b=0) -> None:
        self.r = int(a)
        self.i = int(b)

    def __setr__(self):
        if self.i == 0:
            return f"{self.r}"
        else:
            return f"({self.r},{self.i}i)"
        
    def __repr__(self):
        if self.i == 0:
            return f"{self.r}"
        else:
            return f"({self.r},{self.i}i)"

    def __complex__(self):
        return (complex(self.r,self.i))
    
    def __eq__(self, other: object) -> bool:
        if type(other) is not Gaussian:
            return False
        return self.r == other.r and self.i == other.i
    
    def __ne__(self,other):
        return not self == other
    
    def conjugate(self):  #Returns the conjugate of an Gaussian Integer
        return Gaussian(self.r, -self.i)
    
    def norm(self):
        return self.r**2 + self.i**2
    
    def __pos__(self):
        return self
    
    def add(self, other):
        real = self.r + other.r
        img = self.i + other.i

        return Gaussian(real,img)
    
    def __add__(self, other):

        if type(other) is int:
            return Gaussian(self.r + other, self.i)      
        return self.add(other)
    
    def __radd__(self, other):
        if type(other) is int:
            return Gaussian(self.r + other, self.i)
        return self.add(other)
    
    def __neg__(self):
        return Gaussian(-self.r,-self.i)
    
    def __sub__(self, other):
        if type(other) is int:
            return Gaussian(self.r-other, self.i)
        return self.add(-other)
    
    def floordiv(self,divisor):
        if type(divisor) is int:		
            numerator = (-self if (divisor < 0) else self)		
            denominator = (-divisor if (divisor < 0) else divisor)			
            if denominator == 0:
			
                raise ZeroDivisionError("{0:s} is null!".format(divisor))		
        else:			
            numerator = self*divisor.conjugate()			
            denominator = divisor.norm()	# Recall that denominator >= 0			
            if denominator == 0:				
                raise ZeroDivisionError("{0:s} is null!".format(divisor))		
        candidate_r = numerator.r//denominator
        candidate_i = numerator.i//denominator
		
		# i.e. (candidate_r+1)*denominator-numerator.r < numerator.r-candidate_r*denominator
        if (2*candidate_r+1)*denominator < 2*numerator.r:			
            candidate_r += 1		
		# i.e. (candidate_i+1)*denominator-numerator.i < numerator.i-candidate_i*denominator
        if (2*candidate_i+1)*denominator < 2*numerator.i:		
            candidate_i += 1
		
        return Gaussian(candidate_r,candidate_i)    
    def __floordiv__(self, divisor):
        return self.floordiv(divisor)
    
    def mod(self, divisor):
        return self - divisor * ( self // divisor)
    
    def __mod__(self, divisor):
        return self.mod(divisor)
    
    def __rsub__(self, other):
        if type(other) is int:
            return Gaussian(self.r - other, self.i)
        return self.add(-other)
    
    def multiplication(self, other):
        return Gaussian(self.r * other.r - self.i * other.i, self.r * other.i + other.r * self.i)
    
    def __mul__(self, other):
        return self.multiplication(other)
    
    def normlist(self):
        N = self.norm()
        norms = []
        k = 2
        while N != 1:
            if N % k == 0:
                 N //= k
                 norms.append(k)
            else:
                k += 1
        return norms
    
    def primefactor(self):
        check = Gaussian(1)
        norms = self.normlist()
        pfactor = []
        for p in norms:
            if p==2:
                pfactor.append(Gaussian(1,1))
            elif (p % 4) == 3:
                pfactor.append(Gaussian(p,0))
                norms.remove(p)
            elif (p % 4) == 1:
                prime_candidate = None
                for a in range(1,p):
                    for i in range(p-1,0,-1):
                        if a*a + i*i == p:
                            prime_candidate = Gaussian(max(a,i),min(a,i))
                            break
                    if prime_candidate is not None:
                        break
                
                result = self
                for prime in pfactor:
                    result = result // prime

                if result % prime_candidate == Gaussian(0):
                    pfactor.append(prime_candidate)
                else:
                    pfactor.append(prime_candidate.conjugate())
                
        for prime in pfactor:
            check *= prime
        unit = self // check
        if unit == Gaussian(1,1):
            pfactor[-1] = Gaussian(pfactor[-1].i,pfactor[-1].r)
        elif unit == Gaussian(1,-1):
            pfactor[-1] = pfactor[-1].conjugate()
            pfactor.append(Gaussian(0,-1))
        elif unit == Gaussian(-1,1):
            pfactor.append(Gaussian(0,1))
        elif unit == Gaussian(-1,-1):
            pfactor[-1] = Gaussian(-pfactor[-1].i,-pfactor[-1].r)
        elif unit == Gaussian(1):
            pass
        else:
            pfactor.append(unit)

        return pfactor
        
    def gcd(self, other):
        K = self.primefactor()
        L = other.primefactor()
        M = []
        for value in K:
            if value in L:
                M.append(value)
                L.remove(value)
        total = Gaussian(1)
        for p in M:
            total = total * p
        return total
