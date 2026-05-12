#import "@preview/touying:0.7.3": *
#import themes.university: *

#import "@preview/mannot:0.3.3": *

#set math.mat(delim: "[")
#set math.vec(delim: "[")
#show: university-theme.with(
  aspect-ratio: "16-9",
  config-info(
    title: "Finding Eigenvalues with the QR Algorithm",
    subtitle: "ES.1806 Final Project",
    author: "Chris Liu",
    date: datetime.today(),
    institution: "MIT",
  ),
)

#title-slide()

== Outline

- Naive approach to eigenvalues: $det(A - lambda bb(1)) = 0$
- Similar matrices
- Finding a single eigenvalue
- Orthogonal iteration
- QR algorithm
  - Example
- Improving the QR algorithm

== Similar matrices

$A$ and $B$ are similar if $A = X B X^(-1)$.

#let vv(t) = math.bold(math.upright(t))
#let TT = math.upright("T")

If $A vv(v)_A = lambda_A vv(v)_A$, then

$
   X B X^(-1) vv(v)_A & = lambda_A vv(v)_A \
  => B X^(-1) vv(v)_A & = lambda_A X^(-1) vv(v)_A
$

$lambda_A$ is also an eigenvalue of $B$ with eigenvector $X^(-1) vv(v)_A$, and vice versa.

== Finding a single eigenvalue

We can write a generic vector $vv(x)$ in the eigenbasis: $vv(x) = sum_(i=1)^n c_i vv(v)_i$. If we
repeatedly apply $A$ to $vv(x)$, we get:

$
  A^k vv(x) & = sum_(i=1)^n c_i A^k vv(v)_i \
            & = sum_(i=1)^n c_i lambda_i^k vv(v)_i \
$

Suppose $lambda_1$ is the biggest eigenvalue, then

$
  A^k vv(x) & = lambda_1^k (c_1 vv(v)_1 + markhl(sum_(i=2)^n c_i (lambda_i/lambda_1)^k vv(v)_i, color: #gray))
$

=== Algorithm

+ Start with any $n times 1$ vector $vv(x)$.
+ Multiply $vv(x)$ by $A$.
+ Normalize the vector and go to step 2.

== Orthogonal iteration

Idea: adapt previous algorithm to find multiple eigenvectors at once! \
$=>$ Use orthogonal matrix $P$ to store a basis for the span of the eigenvectors

$
  => P_(k+1)R_(k+1) = A P_k
$

=== Algorithm

+ Choose any $n times n$ matrix $P_0$. Usually we pick $P_0 = bb(1)_n$ for simplicity.
+ Calculate $A P_k$.
+ Decompose the result into its QR factorization $P_(k+1)R_(k+1)$. Go to step 2.

If the algorithm converges, diagonal entries of $R_k$ contain the eigenvalues of $A$:

$
  P_(k+1)R_(k+1) & = A P_k \
         R_(k+1) & = P_(k+1)^(-1) A P_k = P_(k+1)^TT A P_k
$

Convergence gives $R_(k+1) approx R_k$, and $P_(k+1) approx P_k$:

$
  R_k approx P_k^TT A P_k
$

$R_k$ is upper triangular $=>$ Eigenvalues are on the diagonal

== QR Algorithm

Observation: $R_k$ is only similar to $A$ at convergence! What if we found a matrix that is _always_
similar to $A$?

Trick: Let $A_k = P_k^TT A P_k$. Define $Q_(k+1) = P_k^TT P_(k+1)$.

#align(center)[_... recurrence magic omitted (see board if time) ..._]

$
  mark(A_k = Q_(k+1) R_(k+1), color: #red) \
  mark(R_(k+1) Q_(k+1) = A_(k+1), color: #red)
$

#pagebreak()

=== Algorithm

+ Start with $A_0 = A$.
+ Find the QR factorization $A_k = Q_(k+1) R_(k+1)$.
+ Reverse the order to find $A_(k+1) = R_(k+1) Q_(k+1)$. Go to step 2.

$A_k$ is _always_ similar to $A$. At convergence,

$
     & P_k approx P_(k+1) \
  => & Q_(k+1) = P_k^TT P_(k+1) approx bb(1) \
  => & A_k approx R_(k+1)
$

=== Example

$
  A = mat(2, 1; 1, 2)
$

#columns(2)[
  $
    vv(q)_1 & = 1/sqrt(5) vec(2, 1) \
    vv(v)_2 & = vec(1, 2) - 1/5 vec(2, 1) mat(2, 1) vec(1, 2) \
            & = 1/5 vec(-3, 6) \
    vv(q)_2 & = vv(v)_2 / abs(vv(v)_2) = 1/sqrt(5) vec(-1, 2) \
  $
  $
         Q_1 & = 1/sqrt(5) mat(2, -1; 1, 2) \
         R_1 & = Q_1^(-1) A = 1/sqrt(5) mat(5, 4; 0, 3) \
         A_1 & = R_1 Q_1 \
             & = 1/5 mat(14, 3; 3, 6) \
    lambda_1 & approx 14/5 = mark(2.8, color: #red) \
    lambda_2 & approx 6/5 = mark(1.2, color: #red) \
  $
]

#columns(2)[
  $
    Q_2 & = 1/sqrt(205) mat(14, -3; 3, 14) \
    R_2 & = 1/sqrt(205) mat(41, 12; 0, 15) \
  $

  #colbreak()

  $
         A_2 & = R_2 Q_2 = 1/41 mat(122, 9; 9, 42) \
    lambda_1 & approx 122/41 approx mark(2.976, color: #red) \
    lambda_2 & approx 42/41 approx mark(1.024, color: #red) \
  $
]

Python simulation:

#set text(20pt)

#align(center, table(
  columns: 8,
  $k$, ..range(7).map(n => $#n$),
  $lambda_1$,
  $2.00000000$,
  $2.80000000$,
  $2.97560976$,
  $2.99726027$,
  $2.99969521$,
  $2.99996613$,
  $2.99999624$,
  $lambda_2$,
  $2.00000000$,
  $1.20000000$,
  $1.02439024$,
  $1.00273973$,
  $1.00030479$,
  $1.00003387$,
  $1.00000376$,
))

== Improving the QR algorithm (if time)

- Shift-inverting the matrix gives faster convergence:

  $
    (A - s bb(1))^(-1)
  $

- Upper Hessenberg Form reduces computational complexity:
  - Householder reflections
  - Givens rotations

  $
    mat(
      a_11, a_12, a_13, a_14;
      a_21, a_22, a_23, a_24;
      0, a_32, a_33, a_34;
      0, 0, a_43, a_44;
    )
  $

== References

+ Jim Lambers, Advanced Topics in Numerical Linear Algebra,
  https://web.stanford.edu/class/cme335/lecture4
+ David Bindle, Matrix Computations,
  https://www.cs.cornell.edu/~bindel/class/cs6210-f12/notes/lec24.pdf
+ David Bindle, Matrix Computations,
  https://www.cs.cornell.edu/~bindel/class/cs6210-f12/notes/lec25.pdf
+ David Bindle, Matrix Computations,
  https://www.cs.cornell.edu/~bindel/class/cs6210-f12/notes/lec26.pdf
+ Gilbert Strang, Computing Eigenvalues and Singular Values,
  https://ocw.mit.edu/courses/18-065-matrix-methods-in-data-analysis-signal-processing-and-machine-learning-spring-2018/resources/lecture-12-computing-eigenvalues-and-singular-values/
+ Martijn Anthonissen, QR algorithm for computing eigenvalues,
  https://www.youtube.com/watch?v=BfGAmw9qKsM
+ Wikipedia, QR algorithm, https://en.wikipedia.org/wiki/QR_algorithm
+ Wikipedia, Schur decomposition, https://en.wikipedia.org/wiki/Schur_decomposition
