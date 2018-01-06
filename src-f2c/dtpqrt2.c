#line 1 "dtpqrt2.f"
/* dtpqrt2.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

#line 1 "dtpqrt2.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b5 = 1.;
static doublereal c_b17 = 0.;

/* > \brief \b DTPQRT2 computes a QR factorization of a real or complex "triangular-pentagonal" matrix, which 
is composed of a triangular block and a pentagonal block, using the compact WY representation for Q. 
*/

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DTPQRT2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtpqrt2
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtpqrt2
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtpqrt2
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DTPQRT2( M, N, L, A, LDA, B, LDB, T, LDT, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER   INFO, LDA, LDB, LDT, N, M, L */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), T( LDT, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DTPQRT2 computes a QR factorization of a real "triangular-pentagonal" */
/* > matrix C, which is composed of a triangular block A and pentagonal block B, */
/* > using the compact WY representation for Q. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The total number of rows of the matrix B. */
/* >          M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix B, and the order of */
/* >          the triangular matrix A. */
/* >          N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] L */
/* > \verbatim */
/* >          L is INTEGER */
/* >          The number of rows of the upper trapezoidal part of B. */
/* >          MIN(M,N) >= L >= 0.  See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the upper triangular N-by-N matrix A. */
/* >          On exit, the elements on and above the diagonal of the array */
/* >          contain the upper triangular matrix R. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is DOUBLE PRECISION array, dimension (LDB,N) */
/* >          On entry, the pentagonal M-by-N matrix B.  The first M-L rows */
/* >          are rectangular, and the last L rows are upper trapezoidal. */
/* >          On exit, B contains the pentagonal matrix V.  See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B.  LDB >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] T */
/* > \verbatim */
/* >          T is DOUBLE PRECISION array, dimension (LDT,N) */
/* >          The N-by-N upper triangular factor T of the block reflector. */
/* >          See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* >          LDT is INTEGER */
/* >          The leading dimension of the array T.  LDT >= max(1,N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0: successful exit */
/* >          < 0: if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup doubleOTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The input matrix C is a (N+M)-by-N matrix */
/* > */
/* >               C = [ A ] */
/* >                   [ B ] */
/* > */
/* >  where A is an upper triangular N-by-N matrix, and B is M-by-N pentagonal */
/* >  matrix consisting of a (M-L)-by-N rectangular matrix B1 on top of a L-by-N */
/* >  upper trapezoidal matrix B2: */
/* > */
/* >               B = [ B1 ]  <- (M-L)-by-N rectangular */
/* >                   [ B2 ]  <-     L-by-N upper trapezoidal. */
/* > */
/* >  The upper trapezoidal matrix B2 consists of the first L rows of a */
/* >  N-by-N upper triangular matrix, where 0 <= L <= MIN(M,N).  If L=0, */
/* >  B is rectangular M-by-N; if M=L=N, B is upper triangular. */
/* > */
/* >  The matrix W stores the elementary reflectors H(i) in the i-th column */
/* >  below the diagonal (of A) in the (N+M)-by-N input matrix C */
/* > */
/* >               C = [ A ]  <- upper triangular N-by-N */
/* >                   [ B ]  <- M-by-N pentagonal */
/* > */
/* >  so that W can be represented as */
/* > */
/* >               W = [ I ]  <- identity, N-by-N */
/* >                   [ V ]  <- M-by-N, same form as B. */
/* > */
/* >  Thus, all of information needed for W is contained on exit in B, which */
/* >  we call V above.  Note that V has the same form as B; that is, */
/* > */
/* >               V = [ V1 ] <- (M-L)-by-N rectangular */
/* >                   [ V2 ] <-     L-by-N upper trapezoidal. */
/* > */
/* >  The columns of V represent the vectors which define the H(i)'s. */
/* >  The (M+N)-by-(M+N) block reflector H is then given by */
/* > */
/* >               H = I - W * T * W**T */
/* > */
/* >  where W^H is the conjugate transpose of W and T is the upper triangular */
/* >  factor of the block reflector. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dtpqrt2_(integer *m, integer *n, integer *l, doublereal *
	a, integer *lda, doublereal *b, integer *ldb, doublereal *t, integer *
	ldt, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, t_dim1, t_offset, i__1, i__2, 
	    i__3;

    /* Local variables */
    static integer i__, j, p, mp, np;
    extern /* Subroutine */ int dger_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal alpha;
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), dtrmv_(char *, 
	    char *, char *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen, ftnlen, ftnlen), dlarfg_(integer *, doublereal 
	    *, doublereal *, integer *, doublereal *), xerbla_(char *, 
	    integer *, ftnlen);


/*  -- LAPACK computational routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments */

#line 208 "dtpqrt2.f"
    /* Parameter adjustments */
#line 208 "dtpqrt2.f"
    a_dim1 = *lda;
#line 208 "dtpqrt2.f"
    a_offset = 1 + a_dim1;
#line 208 "dtpqrt2.f"
    a -= a_offset;
#line 208 "dtpqrt2.f"
    b_dim1 = *ldb;
#line 208 "dtpqrt2.f"
    b_offset = 1 + b_dim1;
#line 208 "dtpqrt2.f"
    b -= b_offset;
#line 208 "dtpqrt2.f"
    t_dim1 = *ldt;
#line 208 "dtpqrt2.f"
    t_offset = 1 + t_dim1;
#line 208 "dtpqrt2.f"
    t -= t_offset;
#line 208 "dtpqrt2.f"

#line 208 "dtpqrt2.f"
    /* Function Body */
#line 208 "dtpqrt2.f"
    *info = 0;
#line 209 "dtpqrt2.f"
    if (*m < 0) {
#line 210 "dtpqrt2.f"
	*info = -1;
#line 211 "dtpqrt2.f"
    } else if (*n < 0) {
#line 212 "dtpqrt2.f"
	*info = -2;
#line 213 "dtpqrt2.f"
    } else if (*l < 0 || *l > min(*m,*n)) {
#line 214 "dtpqrt2.f"
	*info = -3;
#line 215 "dtpqrt2.f"
    } else if (*lda < max(1,*n)) {
#line 216 "dtpqrt2.f"
	*info = -5;
#line 217 "dtpqrt2.f"
    } else if (*ldb < max(1,*m)) {
#line 218 "dtpqrt2.f"
	*info = -7;
#line 219 "dtpqrt2.f"
    } else if (*ldt < max(1,*n)) {
#line 220 "dtpqrt2.f"
	*info = -9;
#line 221 "dtpqrt2.f"
    }
#line 222 "dtpqrt2.f"
    if (*info != 0) {
#line 223 "dtpqrt2.f"
	i__1 = -(*info);
#line 223 "dtpqrt2.f"
	xerbla_("DTPQRT2", &i__1, (ftnlen)7);
#line 224 "dtpqrt2.f"
	return 0;
#line 225 "dtpqrt2.f"
    }

/*     Quick return if possible */

#line 229 "dtpqrt2.f"
    if (*n == 0 || *m == 0) {
#line 229 "dtpqrt2.f"
	return 0;
#line 229 "dtpqrt2.f"
    }

#line 231 "dtpqrt2.f"
    i__1 = *n;
#line 231 "dtpqrt2.f"
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        Generate elementary reflector H(I) to annihilate B(:,I) */

#line 235 "dtpqrt2.f"
	p = *m - *l + min(*l,i__);
#line 236 "dtpqrt2.f"
	i__2 = p + 1;
#line 236 "dtpqrt2.f"
	dlarfg_(&i__2, &a[i__ + i__ * a_dim1], &b[i__ * b_dim1 + 1], &c__1, &
		t[i__ + t_dim1]);
#line 237 "dtpqrt2.f"
	if (i__ < *n) {

/*           W(1:N-I) := C(I:M,I+1:N)^H * C(I:M,I) [use W = T(:,N)] */

#line 241 "dtpqrt2.f"
	    i__2 = *n - i__;
#line 241 "dtpqrt2.f"
	    for (j = 1; j <= i__2; ++j) {
#line 242 "dtpqrt2.f"
		t[j + *n * t_dim1] = a[i__ + (i__ + j) * a_dim1];
#line 243 "dtpqrt2.f"
	    }
#line 244 "dtpqrt2.f"
	    i__2 = *n - i__;
#line 244 "dtpqrt2.f"
	    dgemv_("T", &p, &i__2, &c_b5, &b[(i__ + 1) * b_dim1 + 1], ldb, &b[
		    i__ * b_dim1 + 1], &c__1, &c_b5, &t[*n * t_dim1 + 1], &
		    c__1, (ftnlen)1);

/*           C(I:M,I+1:N) = C(I:m,I+1:N) + alpha*C(I:M,I)*W(1:N-1)^H */

#line 249 "dtpqrt2.f"
	    alpha = -t[i__ + t_dim1];
#line 250 "dtpqrt2.f"
	    i__2 = *n - i__;
#line 250 "dtpqrt2.f"
	    for (j = 1; j <= i__2; ++j) {
#line 251 "dtpqrt2.f"
		a[i__ + (i__ + j) * a_dim1] += alpha * t[j + *n * t_dim1];
#line 252 "dtpqrt2.f"
	    }
#line 253 "dtpqrt2.f"
	    i__2 = *n - i__;
#line 253 "dtpqrt2.f"
	    dger_(&p, &i__2, &alpha, &b[i__ * b_dim1 + 1], &c__1, &t[*n * 
		    t_dim1 + 1], &c__1, &b[(i__ + 1) * b_dim1 + 1], ldb);
#line 255 "dtpqrt2.f"
	}
#line 256 "dtpqrt2.f"
    }

#line 258 "dtpqrt2.f"
    i__1 = *n;
#line 258 "dtpqrt2.f"
    for (i__ = 2; i__ <= i__1; ++i__) {

/*        T(1:I-1,I) := C(I:M,1:I-1)^H * (alpha * C(I:M,I)) */

#line 262 "dtpqrt2.f"
	alpha = -t[i__ + t_dim1];
#line 264 "dtpqrt2.f"
	i__2 = i__ - 1;
#line 264 "dtpqrt2.f"
	for (j = 1; j <= i__2; ++j) {
#line 265 "dtpqrt2.f"
	    t[j + i__ * t_dim1] = 0.;
#line 266 "dtpqrt2.f"
	}
/* Computing MIN */
#line 267 "dtpqrt2.f"
	i__2 = i__ - 1;
#line 267 "dtpqrt2.f"
	p = min(i__2,*l);
/* Computing MIN */
#line 268 "dtpqrt2.f"
	i__2 = *m - *l + 1;
#line 268 "dtpqrt2.f"
	mp = min(i__2,*m);
/* Computing MIN */
#line 269 "dtpqrt2.f"
	i__2 = p + 1;
#line 269 "dtpqrt2.f"
	np = min(i__2,*n);

/*        Triangular part of B2 */

#line 273 "dtpqrt2.f"
	i__2 = p;
#line 273 "dtpqrt2.f"
	for (j = 1; j <= i__2; ++j) {
#line 274 "dtpqrt2.f"
	    t[j + i__ * t_dim1] = alpha * b[*m - *l + j + i__ * b_dim1];
#line 275 "dtpqrt2.f"
	}
#line 276 "dtpqrt2.f"
	dtrmv_("U", "T", "N", &p, &b[mp + b_dim1], ldb, &t[i__ * t_dim1 + 1], 
		&c__1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*        Rectangular part of B2 */

#line 281 "dtpqrt2.f"
	i__2 = i__ - 1 - p;
#line 281 "dtpqrt2.f"
	dgemv_("T", l, &i__2, &alpha, &b[mp + np * b_dim1], ldb, &b[mp + i__ *
		 b_dim1], &c__1, &c_b17, &t[np + i__ * t_dim1], &c__1, (
		ftnlen)1);

/*        B1 */

#line 286 "dtpqrt2.f"
	i__2 = *m - *l;
#line 286 "dtpqrt2.f"
	i__3 = i__ - 1;
#line 286 "dtpqrt2.f"
	dgemv_("T", &i__2, &i__3, &alpha, &b[b_offset], ldb, &b[i__ * b_dim1 
		+ 1], &c__1, &c_b5, &t[i__ * t_dim1 + 1], &c__1, (ftnlen)1);

/*        T(1:I-1,I) := T(1:I-1,1:I-1) * T(1:I-1,I) */

#line 291 "dtpqrt2.f"
	i__2 = i__ - 1;
#line 291 "dtpqrt2.f"
	dtrmv_("U", "N", "N", &i__2, &t[t_offset], ldt, &t[i__ * t_dim1 + 1], 
		&c__1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*        T(I,I) = tau(I) */

#line 295 "dtpqrt2.f"
	t[i__ + i__ * t_dim1] = t[i__ + t_dim1];
#line 296 "dtpqrt2.f"
	t[i__ + t_dim1] = 0.;
#line 297 "dtpqrt2.f"
    }

/*     End of DTPQRT2 */

#line 302 "dtpqrt2.f"
    return 0;
} /* dtpqrt2_ */

