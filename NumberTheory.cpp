#include <bits/stdc++.h>

using namespace std;

typedef long long int _ll;

#define MOD 1000000007

#define MAX_N 100

_ll facts[MAX_N];
_ll inverseFacts[MAX_N];
_ll phi[MAX_N];

void shieve()
{

  // Shieve for phi starts here
  for (int i = 0; i < MAX_N; i++)
    phi[i] = i;

  phi[1] = 1;

  for (int i = 2; i < MAX_N; i++)
  {
    // i is the divisors of 2*i, 3*i ....
    if (phi[i] == i)
    {
      // Prime case
      phi[i] = i - 1;
      for (_ll j = i + i; j < MAX_N; j += i)
      {
        phi[j] -= phi[j] / i;
      }
    }
  }

  // shieve for phi ends here
}

_ll power(_ll a, _ll p, _ll mod)
{

  _ll result = 1;
  while (p > 0)
  {
    if (p & 1)
    {
      result = (result * a) % mod;
    }

    a = (a * a) % mod;
    p >>= 1;
  }
  return result;
}

void learnFactorial()
{
  facts[0] = 1;
  for (int i = 1; i < MAX_N; i++)
  {
    facts[i] = (i * facts[i - 1]) % MOD;
  }

  // Inverse Factorial
  inverseFacts[MAX_N - 1] = power(facts[MAX_N - 1], MOD - 2, MOD); // Fermat theorem to calculate inverse of a number

  /*
    a! * x === 1 (mod MOD) => find x
  */

  for (int i = MAX_N - 2; i >= 0; i--)
  {
    inverseFacts[i] = (inverseFacts[i + 1] * (i + 1)) % MOD;
  }
}

int gcd(int a, int b)
{
  return b % a ? gcd(b % a, a) : a;
}

void learnGcd()
{

  // GCD (a,b) = greatest common divisor of number a and b
  /*
    Properties :
      GCD(a, GCD(b,c)) = GCD(GCD(a,b), c)
      GCD(a,b) <= min(a,b)
  */

  int a = 12, b = 8;
  cout << "GCD of (" << a << "," << b << ") : " << gcd(a, b) << endl;
}

int extendedEuclidian(int a, int b, int &x, int &y)
{

  if (a == 0)
  {
    //  a*x + b*y = gcd(a,b)
    // b =0, so gcd(a,b) = a, and hence x = 1,
    x = 0;
    y = 1;
    return b;
  }

  int x1, y1, g;
  g = extendedEuclidian(b % a, a, x1, y1);

  x = y1 - (b / a) * x1;
  y = x1;
  return g;
}

void learnExtendedGcd()
{
  /*
    let say two numbers x and y,
    There always exists two number a and b such that ax + by = gcd(x,y)

    recursively we can the coefficient a b and gcd(x,y) = g

    ax + by = g
      if gcd(a,b) = 1 and let say b = n
    ax + ny = 1
    ax === 1 (mod n)
    Here x is nothing but the inverse modular of a with respect to n.

    if n is prime then,
      modInverse of a wrt to n is = power(a, n-2, n);
      Proof:
        a ^ phi(n) === 1 (mod n) if gcd(a,n) = 1
        if n is already prime then gce(a,n) = 1 && phi(n) = n-1;
        a^(n-1) === 1 (mod n)
        a^(-1) === a^(n-2) (mod n) => Multiplying both side by a^-1
                       == power(a,n-2,n)


  */

  int x = 0, y = 0;
  int a = 12, b = 8;

  int g = extendedEuclidian(a, b, x, y);

  cout << "The equation " << a << "x" << x << "+" << b << "x" << y << "=" << g << endl;
}

void learnLcm()
{
  /*
    LCM (a, b) = Least common multiplicative of a and b
    LCM(a,b) >= min(a,b)
    LCM(a,b) <= a*b
    LCM(a,b) = a*b / GCD(a,b)
  */

  int a = 12, b = 8;
  cout << "LCM of (" << a << "," << b << ") : " << (a * b) / gcd(a, b) << endl;
}

vector<int> getAllPrimeFactors(int n)
{

  // Sqrt(N)

  vector<int> ans;

  for (int i = 3; i <= sqrt(n); i += 2)
  {
    if (n % i == 0)
    {
      ans.push_back(i);
      while (n % i == 0)
        n /= i;
    }
  }

  if (n > 1)
  {
    ans.push_back(n);
  }

  return ans;
}

void learnPrimeFactors()
{
  /*
    primeNumbers : A unit, divisible by only 1 and itself
    Composite Numbers: A number which are the product of prime numbers.
    Any composite number can be uniquely represented as product of its prime factors
  */

  int n = 12131;
  cout << "Prime factors of " << n << " : ";
  for (int &i : getAllPrimeFactors(n))
  {
    cout << i << " ";
  }
  cout << endl;
}

int calculateBinomial(int n, int k)
{
  if (k > n)
  {
    return 0;
  }

  if (n == 0 || k == 0)
  {
    return 1;
  }
  // return facts[n]/(facts[n-k] * facts[k]); // => It will have problem when overflow will happen
  return (((facts[n] * inverseFacts[n - k]) % MOD) * inverseFacts[k]) % MOD;
}

int binomialCoefficientWithRecurrence(int n, int k)
{
  if (k > n)
  {
    return 0;
  }
  if (n == 0 || k == 0 || n == k)
  {
    return 1;
  }

  return binomialCoefficientWithRecurrence(n - 1, k - 1) + binomialCoefficientWithRecurrence(n - 1, k);
}

void learnBinomialCoefficient()
{
  /*
    A binomial coefficient C(n, k) is the coefficient of x^k in the expansion of (1+x)^n

    Properties:
      C(n,k) = C(n-1, k-1) + C(n-1, k)
      C(n,0) = C(n,n) = 1

      Binomial triangle
                1                   n = 0
              1   1                 n = 1
            1   2   1               n = 2
          1   3   3   1             n = 3
        1   4   6   4   1           n = 4
      1   5   10  10  5   1         n = 5

    This is same as nCk => n choose k => n!/k!(n-k)!

  */
  int n = 10, k = 4;
  cout << "Using factorial  Binomial coefficient for (" << n << "," << k << ") : " << calculateBinomial(n, k) << endl;
  cout << "Using recurssion Binomial coefficient for (" << n << "," << k << ") : " << binomialCoefficientWithRecurrence(n, k) << endl;
}

void learnCatlanNumber()
{
  /*
    Catalan numbers are a sequence of natural numbers that occurs in many interesting counting problems like following.
      Count the number of expressions containing n pairs of parentheses which are correctly matched. For n = 3, possible expressions are ((())), ()(()), ()()(), (())(), (()()).
      Count the number of possible Binary Search Trees with n keys (See this)
      Count the number of full binary trees (A rooted binary tree is full if every vertex has either two children or no children) with n+1 leaves.
      Given a number n, return the number of ways you can draw n chords in a circle with 2 x n points such that no 2 chords intersect.

      Catalan numbers satisfy the following recursive formula.
        C0 = 1
        Cn+1 = Sum (Ci.Cn-i) for all i >= 0 && i <= n
      First few Catlan Numbers
        1, 1, 2, 5, 14, 42, 132

      Formula to calculate nth catlan number
        Cn = (1/n+1)*C(2n, n)

  */
  int n = 10;
  int catlanNumbers[n];
  catlanNumbers[0] = 1;
  cout << "First " << n << " Catlan Numbers are 1 ";
  for (int i = 1; i < n; i++)
  {
    int sum = 0;
    for (int j = 0; j < i; j++)
    {
      sum += catlanNumbers[j] * catlanNumbers[i - j - 1];
    }
    catlanNumbers[i] = sum;
    cout << sum << " ";
  }

  cout << endl;
}

void learnEuclidLemma()
{
  /*
    Euclid’s lemma states that if a prime p divides the product of two numbers (x*y), it must divide at least one of those numbers.
    For example x = 15, y = 6 and p = 5. p divides the product 15*6, it also divides 15.
    The idea is simple, since p is prime, it cannot be factorized. So it must either be completely present in x or in y.
  */
}

void learnEulerTotientFunction()
{
  /*
    Euler’s Totient function Φ (n) for an input n is the count of numbers in {1, 2, 3, …, n} that are relatively prime to n, i.e., the numbers whose GCD (Greatest Common Divisor) with n is 1.
    phi(N) = N*(1-1/p1)*(1-1/p2)......*(1-p/pn)

    Properties:
      phi(a*b) = phi(a) * phi(b) when gcd(a,b) = 1
      phi(p) = p-1 where p is prime (self understood)
      phi(p*q) =  phi(p) * phi(q) =  (p-1) * (q-1) where p and q is prime, p != q
      phi(p^k) = p^k - p^(k-1)
      sum of phi of all divisors of a number n is also n.
          E.g.
            n = 6
            divisors of 6 = 1, 2, 3, 6
                          = phi(1) + phi(2) + phi(3) + phi(6)
                          = 1 + 1 + 2 + 2 = 6

      If  gcd(a, n) == 1,
        then a ^ (phi(n)) === 1 (mod n)
          extension:
            in n is prime => a ^ (n-1) === 1 (mod n)
        E.g.
          a = 2 , n = 5
          2^phi(5) = 16 === 1 (mod 5)
          a = 4 , b = 9
          4 ^(phi(9)) = 4^ 6 = 4096 % 9 === 1 (mod 9)

        => RSA algorithm uses this theorem
  */

  int n = 9;
  int temp = n;
  int res = n;

  for (int i = 2; i <= sqrt(n); i++)
  {
    if (n % i == 0)
    {
      res -= res / i;
      while (n % i == 0)
      {
        n /= i;
      }
    }
  }

  if (n > 1)
  {
    res -= res / n;
  }

  cout << "From sqrt(N) Phi (" << temp << ")=" << res << endl;

  // Check the shieve function to calcaulate phi(n) for all 1 to N

  cout << "From shieve Phi (" << temp << ")=" << phi[temp] << endl;
}

void learnMultiplicativeOrder()
{
  /*
    Let say two number A, N such that gcd(A,N) = 1 , then least k (> 0) which satisfies
      A^k === (1 mod N)
        K is called multiplicative order of A with respect to N
  */

  int A = 4, N = 7;

  int K = 1;
  int res = 1;
  while (K < N)
  {
    res = res * A;
    if (res % N == 1)
    {
      break;
    }
    K++;
  }

  cout << "Multiplicative order of " << A << " with respect to " << N << " is  " << K << endl;
}

void learnNCR()
{

  /*
    Formula:
      nCr = C(n,r) = C(n-1, r-1) + C(n-1, r)
      c(n,0) = C(n,n) = 1
      C(n,r) = n!/r!(n-r)!
      C(n,r)=C(n,r-1)* (n-r+1)/r
        How this formula works
          (A + B)^n = A^nB^0  + nA^n-1B^1 + n*(n-1)/2 A^n-2B^2  ===> Term i, Coefficient of previous term * (power of A in ith term)/total number of terms => C(n,r-1)*(n-(r-1))/r = C(n,r-1)*(n-r+1)/r
    */

  int n = 5, r = 3;
  r = min(n - r, r); // C(n,r) = C(n, n-r)

  int prev = 1;

  // Base case

  for (int i = 1; i <= r; i++)
  {
    prev = prev * (n - i + 1) / i;
  }

  cout << "Value of C(" << n << "," << r << ") : " << prev << endl;
}

int nCrSimple(int n, int r, int p)
{

  r = min(n - r, r); // C(n,r) = C(n, n-r)

  int prev = 1;
  // Base case

  for (int i = 1; i <= r; i++)
  {
    prev = prev * (n - i + 1) / i;
    prev %= p;
  }

  return prev;
}

int nCrByLucas(int n, int r, int p)
{
  if (r == 0)
  {
    return 1;
  }
  if (n == r)
  {
    return 1;
  }

  // convert n in base p, and r in base p
  /*
    Take n modulo p, n0, Divide n with p = n
    Take n modulo p, n1, Divide n with p = n
    till n or r == 0, since n > r, for sure r will become 0 first

  */
  return nCrByLucas(n / p, r / p, p) * nCrSimple(n % p, r % p, p);
}

void learnLucasTheorem()
{
  /*
    Use case,
      Find number C(n,r)%p

      Lucas theorem :
        For non-negative n, r and a prime p, following holds
          C(n,r)%p = Product {for ( i = 0 to k), C(ni, ri)} % p
          where
            n = nk * p^k + n(k-1) * p ^ (k-1) + ................................ + n1 * p + n0
            r = rk * p^k + ..................................................... + r1 * p + r0
        Basically says that,
          To calculate C(n,r)%p can be computed as prodcut of C(ni,ri)
            where each ni, ri are the digits when n is represented in base p, and r is represented in base p

      For e.g.
        C(12,4) % 5 => (12*11*10*9)/(4*3*2*1) = 495 % 5 = 0
          12 = 22 => In base 5
          4  = 04 => in base 5

        C(12,4) = C(2,0) * C(2,4) % 5
                = 1 * 0 % 5 = 0
  */
  int n = 12, r = 4, p = 13;
  cout << "By Lucas theorem : C(" << n << "," << r << ") : " << nCrByLucas(n, r, p) << endl;
  cout << "By simple        : C(" << n << "," << r << ") : " << nCrSimple(n, r, p) << endl;
}

void learnFermatTheorem()
{
  /*
    Given three numbers n, r and p, compute the value of nCr mod p. Here p is a prime number greater than n.
    Fermat little theorem says
      if p is a prime number, then for any integer a, the number a^p – a is an integer multiple of p. In the notation of modular arithmetic, this is expressed as:
        a^p = a (mod p)
        Proof comes from EulerTotientFunction
          a ^ phi(p) === 1 (mod p)
          a ^(p-1) === 1 (mod p)
          a^p === a (mod p)
          a^p - a === 0 (mod p)
          This says that p must be a factor of a^p - a.
    nCr can be calculated with this. n! * modInverse(r!) * modInverse (n-r!)
  */
}

void learnChineaseRemainderTheorem()
{
  /*
    Used to solve the system of linear congrunece
      x === r1 ( mod n1)
      x === r2 ( mod n2)
      x === r3 (mod n3)
      x === r4 ( mod n4)
      .
      .
      .
      x === rk ( mod nk)
        where n1, n2 ........ nk are pairwise coprime. x always exists.
      Find x ?

    Chinease Remainder theorem :
      let n = Product{i = 1 to k , ni}
        nn[i] = n/n[i]
        rem[i] = [r1, r2.......rk]
        inv[i] = modInverse of nn[i] with respect to n[i]

      then x = Sum {i = 1 to k , rem[i]*nn[i]*inv[i] } % n
  */

  int num[] = {3, 4, 5}, rem[] = {2, 3, 1};
  int N = sizeof(num)/sizeof(int);

  int n = 1;
  for ( int i = 0 ; i < N ; i++){
    n *= num[i];
  }
  int x = 0;
  for ( int i = 0 ; i < N ; i++ ){
    int nn = n/num[i];
    int x1 = 0 , y1 = 0;
    int g = extendedEuclidian(nn, num[i], x1, y1);
    x += rem[i] * nn * x1 ;
    x %= n;
  }

  cout << (x + n)%n << endl ;

}

void squareRoot(int n , int p){

  if ( p % 4 != 3){
    cout << "Invlaid Input\n";
    return;
  }

  // x = +n^(p+1)/4
  n %= p;
  int x = power(n , (p+1)/4 , p);
  if ( (x*x) % p == n ){
    cout << "Square Root of " << n << " mod " << p << " is : " << x << endl ;
    return;
  }

  // Try x = -n^(p+1)/4
  x = p-x ;
  if ((x*x) % p == n){
    cout << "Square Root of " << n << " mod " << p << " is : " << x << endl ;
    return;
  }

  cout << "Square root doesn't exist \n";

}

void squareRootUnderModulo(){

    /*
      Given a number ‘n’ and a prime ‘p’, find square root of n under modulo p if it exists. It may be given that p is in the form for 4*i + 3 (OR p % 4 = 3) where i is an integer. Examples of such primes are 7, 11, 19, 23, 31, … etc,
      Input:  n = 2, p = 7
      Output: 3 or 4
      3 and 4 both are square roots of 2 under modulo
      7 because (3*3) % 7 = 2 and (4*4) % 7 = 2

      Find x such that
        x^2 === n mod p

      Brute Force:
        take number i from 1 to p-1, check for all i*i is n mod p or not.

      if p is in the form of 4n+3,

        If square root exists, then
          n^(p-1)/2 === 1 (mod p)
          n^(p+1)/2 === n (mod p) => Multiplying both side by n

        Let x be the modulo square root
            (x * x) ≡ n mod p
            (x * x) ≡ n^(p+1)/2  [Using (1) given above]
            (x * x) ≡ n^(2i + 2) [Replacing n = 4*i + 3]
            x ≡ ±n^(i + 1)  [Taking Square root of both sides]
            x ≡ ±n^(p + 1)/4 [Putting 4*i + 3 = p or i = (p-3)/4]
    */

    squareRoot(2,7);

}

int main()
{

  ios::sync_with_stdio(false);
  cin.tie(NULL);
  cout.tie(NULL);

  shieve();
  learnFactorial();
  learnGcd();
  learnExtendedGcd();
  learnLcm();
  learnPrimeFactors();
  learnBinomialCoefficient();
  learnCatlanNumber();
  learnEuclidLemma();

  learnEulerTotientFunction();
  learnMultiplicativeOrder();

  learnNCR();
  learnLucasTheorem();
  learnFermatTheorem();

  learnChineaseRemainderTheorem();

  squareRootUnderModulo();
}