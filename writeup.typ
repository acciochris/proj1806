#import "@preview/academic-alt:0.1.0": *
#import "@preview/lemmify:0.1.8": *

#set page(paper: "us-letter")
#set math.equation(numbering: "(1)")
#set math.mat(delim: "[")
#set math.vec(delim: "[")
#let (
  theorem,
  lemma,
  corollary,
  remark,
  proposition,
  example,
  proof,
  rules: thm-rules,
) = default-theorems("thm-group", lang: "en", thm-numbering: thm-numbering-linear)
#show: thm-rules
#show: university-assignment.with(
  title: "Finding Eigenvalues with\nthe QR Algorithm",
  subtitle: "ES.1806 Final Project",
  author: "Chris Liu",
  details: (
    course: "ES.1806 Linear Algebra",
    instructor: "Arthur Parzygnat",
  ),
)

= Introduction

Eigenvalues and eigenvectors are an immensely useful aspect of linear algebra in that they have
wide-ranging applications such as solving differential equations and finding high powers of matrices
with little computation. However, the naive approach of finding eigenvalues, namely, solving the
equation $det(A - lambda bb(1)) = 0$ is computationally expensive and analytically impossible for
polynomials of degree 5 or higher. In this project, we investigate an alternative method of
computing eigenvalues and eigenvectors numerically: the QR algorithm. As the name suggests, the
algorithms involves repeatedly finding the QR decomposition of a matrix.

Before we describe the QR algorithm, we need to explain two simple concepts, both of which are
extensively used in the derivation of the QR algorithm.

== Similar matrices

Two square $n times n$ matrices $A$ and $B$ are similar if there exists an invertible matrix $X$
such that $A = X B X^(-1)$.

#let vv(t) = math.bold(math.upright(t))
#lemma[If $A$ and $B$ are similar, they have the same set of eigenvalues.] <lemma1>
#proof[
  Suppose $lambda_A$ is an eigenvalue of $A$, then by definition there exists $vv(v)_A$ such that
  $A vv(v)_A = lambda_A vv(v)_A$. Plugging in $A = X B X^(-1)$, we get

  $
     X B X^(-1) vv(v)_A & = lambda_A vv(v)_A \
    => B X^(-1) vv(v)_A & = lambda_A X^(-1) vv(v)_A
  $

  Therefore, $lambda_A$ is an eigenvalue of $B$ with eigenvector $X^(-1) vv(v)_A$. If we replace $X$
  with its inverse, we can make an identical argument to prove that every eigenvalue of $B$ is also
  an eigenvalue of $A$.
]

== Finding a single eigenvalue <eigen1>

While there is no simple method of finding all eigenvalues of an arbitrarily large matrix. There is
a simple algorithm to find the eigenvalue with the largest magnitude. For the sake of simplicity, we
assume the $n times n$ matrix $A$ is has $n$ _real, positive, distinct_ eigenvalues. We also assume
the eigenvalues $lambda_1, lambda_2, ..., lambda_n$ are arranged in strictly decreasing order, and
that their corresponding eigenvectors are $vv(v)_1, vv(v)_1, ..., vv(v)_n$. Suppose the eigenvalue
decomposition $A = X Lambda X^(-1)$ also follows this order.

Suppose we have a generic $n times 1$ vector $vv(x)$. We can rewrite the quantity $A^k vv(x)$, i.e.
$A$ applied to $vv(x)$ $k$ times as follows:

$
  A^k vv(x) = (X Lambda X^(-1))^k vv(x) = X Lambda^k X^(-1) vv(x)
$

What does this expression mean? $X^(-1) vv(x)$ is the coordinates of $vv(x)$ in the eigenbasis, so
$X Lambda^k X^(-1) vv(x)$ would be the linear combination of the eigenvectors of $A$ where the
coefficient of each $vv(v)_i$ is given by the product of $lambda_i^k$ and the $i$th coordinate of
$vv(x)$ in the eigenbasis. To put it more clearly, suppose

$
  vv(x) = sum_(i=1)^n c_i vv(v)_i = X vec(c_1, dots.v, c_n),
$

then

$
  A^k vv(x) = X Lambda^k X^(-1) X vec(c_1, dots.v, c_n)
  = X mat(lambda_1^k, "", ""; "", dots.down, ""; "", "", lambda_n^k) vec(c_1, dots.v, c_n)
  = X vec(lambda_1^k c_1, dots.v, lambda_n^k c_n)
  = sum_(i=1)^n lambda_i^k c_i vv(v)_i
$

We can rescale the resulting vector by $1/lambda_1^k$, where $lambda_1$ is the largest eigenvalue:

$
  1/lambda_1^k A^k vv(x) = c_1 vv(v)_1 + sum_(i=2)^n (lambda_i/lambda_1)^k c_i vv(v)_i
$

because $lambda_i / lambda_1 < 1$ for all $i >= 2$, the second term disappears as we take the limit
$k -> infinity$. The only remaining term is $c_1 vv(v)_1$, and we have now found an eigenvector of
$A$!

To summarize, the following algorithm finds the eigenvector with the largest eigenvalue for an
$n times n$ matrix $A$ (with several limitations described at the beginning of this section):

+ Start with any $n times 1$ vector $vv(x)$.
+ Multiply $vv(x)$ by $A$.
+ Normalize the vector and go to step 2. (We technically don't need the normalization step, but it
  is difficult to determine convergence if the magnitude of the vector keeps changing.)

After an adequate number of iterations, the resulting vector will be a decent approximation of
$vv(v)_1$ and we can find the corresponding eigenvalue by using $A vv(v)_1 = lambda_1 vv(v)_1$.


= The QR algorithm

We motivate the QR algorithm by first looking at orthogonal iteration, which is an extension of the
algorithm in @eigen1. We then discuss how the canonical QR algorithm is merely a clever manipulation
of the orthogonal iteration algorithm using recurrence relationships.

== Orthogonal iteration

#let span = math.op("span")
#let CC = math.op(math.bold("C"))
#let TT = math.upright("T")

The basic idea of orthogonal iteration is that instead of finding a single eigenvector of $A$, we
find the subspace spanned by the eigenvectors, i.e. $CC(X)$. If we maintain the same assumptions as
@eigen1 (real, positive, distinct eigenvalues), $A$ is diagonalizable, and we have $CC(X) = CC(A)$.
In other words, the eigenvectors span the column space of $A$.

To actually implement the algorithm, we need to maintain a matrix whose columns serve as the basis
of the vector space we're trying to find. Specifically, we maintain an orthogonal matrix $P_k$ (it
is *not* called $Q_k$ for a reason that will become evident in @recur). The subscript $k$ refers to
the $k$th iteration of the algorithm.

The algorithm is as follows:

+ Choose any $n times n$ matrix $P_0$. Usually we pick $P_0 = bb(1)_n$ for the purpose of
  simplicity.
+ Calculate $A P_k$.
+ Decompose the result into its QR factorization $P_(k+1)R_(k+1)$. (Once again, I apologize for the
  fact that the orthogonal matrix is not called $Q$.) Go to step 2.

Mathematically, the following recurrence relationship is maintained:

$
  P_(k+1)R_(k+1) = A P_k
$ <eq1>

Notice that the same principles apply as in section @eigen1. We repeatedly apply $A$ to an entity
and then rescale it in some way after each iteration, in this case using the QR factorization. In
fact, if we only consider the first column $vv(p)$ of matrix $P$ on both sides of @eq1, we have:

$
  vv(p)_(k+1)r_(k+1) = A vv(p)_k
$ <eq2>

($r_(k+1)$ is a scalar because the matrix $R_(k+1)$ is upper triangular.) If we squint our eyes at
@eq2, we will notice that it is precisely the same algorithm as @eigen1: repeatedly apply $A$ to a
vector, and then normalize it!

A partial proof of correctness of the orthogonal iteration algorithm is given as follows:

#theorem[
  If the orthogonal iteration algorithm converges, the diagonal entries of $R_k$ contain the
  eigenvalues of $A$.
] <thm1>
#proof[
  Because $P_(k+1)$ is orthogonal, we can rewrite @eq1 as:

  $
    R_(k+1) = P_(k+1)^(-1) A P_k = P_(k+1)^TT A P_k
  $

  Convergence gives $R_(k+1) approx R_k$, and $P_(k+1) approx P_k$, which gives

  $
    R_k approx P_k^TT A P_k
  $

  In other words, $R_k$ is similar to $A$. @lemma1 shows that $R_k$ has the same eigenvalues as $A$.
  Since $R_k$ is upper triangular, its eigenvalues are simply its diagonal entries, because the
  determinant of an upper triangular matrix is equal to the product of its diagonals.
]

Proving the convergence of the orthogonal iteration algorithm is beyond the scope of this project.

== Recurrence magic <recur>

Notice that in the orthogonal iteration algorithm, the orthogonal matrix $P_k$ is what gets passed
around between two consecutive iterations. However, $P_k$ does not necessarily contain any useful
information other than the fact that it is orthogonal and that it spans the column space of $A$.
Specifically, it does not contain the eigenvectors of $A$.

In order to find the eigenvalues of $A$, we need to extract them from $R_k$, which is _only_ equal
to $P_k^TT A P_k$ _at convergence_. The canonical version of the QR algorithm mitigates this
inconvenience by translating the orthogonal iteration algorithm into a recurrence relationship that
produces a matrix which is _always_ similar to $A$. As a result, we can keep executing this
algorithm and directly read off the answer once the matrix has converged.

Define the matrix

$
  A_k = P_k^TT A P_k
$

(Note that this is true by definition regardless of whether $P_k$ has converged.)

We can plug in @eq1 to get:

$
  A_k = P_k^TT (A P_k) = P_k^TT P_(k+1) R_(k+1)
$

Let us define the orthogonal matrix $Q_(k+1) = P_k^TT P_(k+1)$ using the fact that the product of
two orthogonal matrices is also orthogonal ($(A B)^(-1) = B^(-1) A^(-1) = B^TT A^TT = (A B)^TT$).
Therefore, $A_k$ itself can now be written as a QR factorization:

$
  A_k = Q_(k+1) R_(k+1)
$

Now consider $A_(k+1)$ and attempt to form a recurrence relation:

$
  A_(k+1) & = P_(k+1)^TT A P_(k+1) \
          & = P_(k+1)^TT P_k A_k P_k^TT P_(k+1) \
          & = Q_(k+1)^TT A_k Q_(k+1) \
          & = Q_(k+1)^TT Q_(k+1) R_(k+1) Q_(k+1) \
          & = R_(k+1) Q_(k+1)
$

Putting everything together, we have:

$
  A_k = Q_(k+1) R_(k+1) \
  R_(k+1) Q_(k+1) = A_(k+1)
$ <eq3>

@eq3 is the essence of the QR algorithm:

+ Start with $A_0 = A$.
+ Find the QR factorization $A_k = Q_(k+1) R_(k+1)$.
+ Reverse the order to find $A_(k+1) = R_(k+1) Q_(k+1)$. Go to step 2.

If the algorithm converges, we have $Q_(k+1) = P_k^TT P_(k+1) approx P_k^TT P_k = bb(1)_n$.
Therefore $A_k approx R_k$, and we can read off the diagonal entries of $A_k$ to find the
eigenvalues by @thm1.

== Example

Consider

$
                        A & = mat(2, 1; 1, 2) \
  det(A - lambda bb(1)_2) & = (2 - lambda)^2 - 1 = 0 \
              => lambda_1 & = 3, lambda_2 = 1
$

Now let's run the QR algorithm. Let $A_0 = A$. Gram-Schmidt gives:

#columns(2)[
  $
    vv(q)_1 & = 1/sqrt(5) vec(2, 1) \
    vv(v)_2 & = vec(1, 2) - 1/5 vec(2, 1) mat(2, 1) vec(1, 2) \
            & = 1/5 vec(-3, 6) \
    vv(q)_2 & = vv(v)_2 / abs(vv(v)_2) = 1/sqrt(5) vec(1, -2) \
  $
  $
         Q_1 & = 1/sqrt(5) mat(2, 1; 1, -2) \
         R_1 & = Q_1^(-1) A = 1/sqrt(5) mat(5, 4; 0, -3) \
         A_1 & = R_1 Q_1 \
             & = 1/sqrt(5) mat(5, 4; 0, -3) 1/sqrt(5) mat(2, 1; 1, -2) \
             & = 1/5 mat(14, -3; -3, 6) \
    lambda_1 & approx 14/5 = 2.8 \
    lambda_2 & approx 6/5 = 1.2 \
  $
]

A second iteration gives:

$
       Q_2 & = 1/sqrt(205) mat(14, 3; -3, 14) \
       R_2 & = 1/sqrt(205) mat(41, -12; 0, 15) \
       A_2 & = R_2 Q_2 \
           & = 1/41 mat(122, -9; -9, 42) \
  lambda_1 & approx 122/41 approx 2.976 \
  lambda_2 & approx 42/41 approx 1.024 \
$

Computer simulation using Python gives the following eigenvalues for further iterations:

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

As we can see, the two eigenvalues converge fairly quickly to their precise values $3$ and $1$.

= Improving the QR algorithm

= References

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
