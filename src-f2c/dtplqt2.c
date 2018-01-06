#line 1 "dtplqt2.f"
/* dtplqt2.f -- translated by f2c (version 20100827).
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

#line 1 "dtplqt2.f"
/* Table of constant values */

static doublereal c_b4 = 1.;
static doublereal c_b10 = 0.;

/* > \brief \b DTPLQT2 computes a LQ factorization of a real or complex "triangular-pentagonal" matrix, which 
is composed of a triangular block and a pentagonal block, using the compact WY representation for Q. 
*/

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DTPLQT2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtplqt2
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtplqt2
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtplqt2
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DTPLQT2( M, N, L, A, LDA, B, LDB, T, LDT, INFO ) */

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
/* > DTPLQT2 computes a LQ a factorization of a real "triangular-pentagonal" */
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
/* >          The number of rows of the lower trapezoidal part of B. */
/* >          MIN(M,N) >= L >= 0.  See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the lower triangular M-by-M matrix A. */
/* >          On exit, the elements on and below the diagonal of the array */
/* >          contain the lower triangular matrix L. */
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
/* >          On entry, the pentagonal M-by-N matrix B.  The first N-L columns */
/* >          are rectangular, and the last L columns are lower trapezoidal. */
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
/* >          T is DOUBLE PRECISION array, dimension (LDT,M) */
/* >          The N-by-N upper triangular factor T of the block reflector. */
/* >          See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* >          LDT is INTEGER */
/* >          The leading dimension of the array T.  LDT >= max(1,M) */
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

/* > \date December 2016 */

/* > \ingroup doubleOTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The input matrix C is a M-by-(M+N) matrix */
/* > */
/* >               C = [ A ][ B ] */
/* > */
/* > */
/* >  where A is an lower triangular N-by-N matrix, and B is M-by-N pentagonal */
/* >  matrix consisting of a M-by-(N-L) rectangular matrix B1 left of a M-by-L */
/* >  upper trapezoidal matrix B2: */
/* > */
/* >               B = [ B1 ][ B2 ] */
/* >                   [ B1 ]  <-     M-by-(N-L) rectangular */
/* >                   [ B2 ]  <-     M-by-L lower trapezoidal. */
/* > */
/* >  The lower trapezoidal matrix B2 consists of the first L columns of a */
/* >  N-by-N lower triangular matrix, where 0 <= L <= MIN(M,N).  If L=0, */
/* >  B is rectangular M-by-N; if M=L=N, B is lower triangular. */
/* > */
/* >  The matrix W stores the elementary reflectors H(i) in the i-th row */
/* >  above the diagonal (of A) in the M-by-(M+N) input matrix C */
/* > */
/* >               C = [ A ][ B ] */
/* >                   [ A ]  <- lower triangular N-by-N */
/* >                   [ B ]  <- M-by-N pentagonal */
/* > */
/* >  so that W can be represented as */
/* > */
/* >               W = [ I ][ V ] */
/* >                   [ I ]  <- identity, N-by-N */
/* >                   [ V ]  <- M-by-N, same form as B. */
/* > */
/* >  Thus, all of information needed for W is contained on exit in B, which */
/* >  we call V above.  Note that V has the same form as B; that is, */
/* > */
/* >               W = [ V1 ][ V2 ] */
/* >                   [ V1 ] <-     M-by-(N-L) rectangular */
/* >                   [ V2 ] <-     M-by-L lower trapezoidal. */
/* > */
/* >  The rows of V represent the vectors which define the H(i)'s. */
/* >  The (M+N)-by-(M+N) block reflector H is then given by */
/* > */
/* >               H = I - W**T * T * W */
/* > */
/* >  where W^H is the conjugate transpose of W and T is the upper triangular */
/* >  factor of the block reflector. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dtplqt2_(integer *m, integer *n, integer *l, doublereal *
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


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

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

#line 212 "dtplqt2.f"
    /* Parameter adjustments */
#line 212 "dtplqt2.f"
    a_dim1 = *lda;
#line 212 "dtplqt2.f"
    a_offset = 1 + a_dim1;
#line 212 "dtplqt2.f"
    a -= a_offset;
#line 212 "dtplqt2.f"
    b_dim1 = *ldb;
#line 212 "dtplqt2.f"
    b_offset = 1 + b_dim1;
#line 212 "dtplqt2.f"
    b -= b_offset;
#line 212 "dtplqt2.f"
    t_dim1 = *ldt;
#line 212 "dtplqt2.f"
    t_offset = 1 + t_dim1;
#line 212 "dtplqt2.f"
    t -= t_offset;
#line 212 "dtplqt2.f"

#line 212 "dtplqt2.f"
    /* Function Body */
#line 212 "dtplqt2.f"
    *info = 0;
#line 213 "dtplqt2.f"
    if (*m < 0) {
#line 214 "dtplqt2.f"
	*info = -1;
#line 215 "dtplqt2.f"
    } else if (*n < 0) {
#line 216 "dtplqt2.f"
	*info = -2;
#line 217 "dtplqt2.f"
    } else if (*l < 0 || *l > min(*m,*n)) {
#line 218 "dtplqt2.f"
	*info = -3;
#line 219 "dtplqt2.f"
    } else if (*lda < max(1,*m)) {
#line 220 "dtplqt2.f"
	*info = -5;
#line 221 "dtplqt2.f"
    } else if (*ldb < max(1,*m)) {
#line 222 "dtplqt2.f"
	*info = -7;
#line 223 "dtplqt2.f"
    } else if (*ldt < max(1,*m)) {
#line 224 "dtplqt2.f"
	*info = -9;
#line 225 "dtplqt2.f"
    }
#line 226 "dtplqt2.f"
    if (*info != 0) {
#line 227 "dtplqt2.f"
	i__1 = -(*info);
#line 227 "dtplqt2.f"
	xerbla_("DTPLQT2", &i__1, (ftnlen)7);
#line 228 "dtplqt2.f"
	return 0;
#line 229 "dtplqt2.f"
    }

/*     Quick return if possible */

#line 233 "dtplqt2.f"
    if (*n == 0 || *m == 0) {
#line 233 "dtplqt2.f"
	return 0;
#line 233 "dtplqt2.f"
    }

#line 235 "dtplqt2.f"
    i__1 = *m;
#line 235 "dtplqt2.f"
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        Generate elementary reflector H(I) to annihilate B(I,:) */

#line 239 "dtplqt2.f"
	p = *n - *l + min(*l,i__);
#line 240 "dtplqt2.f"
	i__2 = p + 1;
#line 240 "dtplqt2.f"
	dlarfg_(&i__2, &a[i__ + i__ * a_dim1], &b[i__ + b_dim1], ldb, &t[i__ *
		 t_dim1 + 1]);
#line 241 "dtplqt2.f"
	if (i__ < *m) {

/*           W(M-I:1) := C(I+1:M,I:N) * C(I,I:N) [use W = T(M,:)] */

#line 245 "dtplqt2.f"
	    i__2 = *m - i__;
#line 245 "dtplqt2.f"
	    for (j = 1; j <= i__2; ++j) {
#line 246 "dtplqt2.f"
		t[*m + j * t_dim1] = a[i__ + j + i__ * a_dim1];
#line 247 "dtplqt2.f"
	    }
#line 248 "dtplqt2.f"
	    i__2 = *m - i__;
#line 248 "dtplqt2.f"
	    dgemv_("N", &i__2, &p, &c_b4, &b[i__ + 1 + b_dim1], ldb, &b[i__ + 
		    b_dim1], ldb, &c_b4, &t[*m + t_dim1], ldt, (ftnlen)1);

/*           C(I+1:M,I:N) = C(I+1:M,I:N) + alpha * C(I,I:N)*W(M-1:1)^H */

#line 253 "dtplqt2.f"
	    alpha = -t[i__ * t_dim1 + 1];
#line 254 "dtplqt2.f"
	    i__2 = *m - i__;
#line 254 "dtplqt2.f"
	    for (j = 1; j <= i__2; ++j) {
#line 255 "dtplqt2.f"
		a[i__ + j + i__ * a_dim1] += alpha * t[*m + j * t_dim1];
#line 256 "dtplqt2.f"
	    }
#line 257 "dtplqt2.f"
	    i__2 = *m - i__;
#line 257 "dtplqt2.f"
	    dger_(&i__2, &p, &alpha, &t[*m + t_dim1], ldt, &b[i__ + b_dim1], 
		    ldb, &b[i__ + 1 + b_dim1], ldb);
#line 259 "dtplqt2.f"
	}
#line 260 "dtplqt2.f"
    }

#line 262 "dtplqt2.f"
    i__1 = *m;
#line 262 "dtplqt2.f"
    for (i__ = 2; i__ <= i__1; ++i__) {

/*        T(I,1:I-1) := C(I:I-1,1:N) * (alpha * C(I,I:N)^H) */

#line 266 "dtplqt2.f"
	alpha = -t[i__ * t_dim1 + 1];
#line 268 "dtplqt2.f"
	i__2 = i__ - 1;
#line 268 "dtplqt2.f"
	for (j = 1; j <= i__2; ++j) {
#line 269 "dtplqt2.f"
	    t[i__ + j * t_dim1] = 0.;
#line 270 "dtplqt2.f"
	}
/* Computing MIN */
#line 271 "dtplqt2.f"
	i__2 = i__ - 1;
#line 271 "dtplqt2.f"
	p = min(i__2,*l);
/* Computing MIN */
#line 272 "dtplqt2.f"
	i__2 = *n - *l + 1;
#line 272 "dtplqt2.f"
	np = min(i__2,*n);
/* Computing MIN */
#line 273 "dtplqt2.f"
	i__2 = p + 1;
#line 273 "dtplqt2.f"
	mp = min(i__2,*m);

/*        Triangular part of B2 */

#line 277 "dtplqt2.f"
	i__2 = p;
#line 277 "dtplqt2.f"
	for (j = 1; j <= i__2; ++j) {
#line 278 "dtplqt2.f"
	    t[i__ + j * t_dim1] = alpha * b[i__ + (*n - *l + j) * b_dim1];
#line 279 "dtplqt2.f"
	}
#line 280 "dtplqt2.f"
	dtrmv_("L", "N", "N", &p, &b[np * b_dim1 + 1], ldb, &t[i__ + t_dim1], 
		ldt, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*        Rectangular part of B2 */

#line 285 "dtplqt2.f"
	i__2 = i__ - 1 - p;
#line 285 "dtplqt2.f"
	dgemv_("N", &i__2, l, &alpha, &b[mp + np * b_dim1], ldb, &b[i__ + np *
		 b_dim1], ldb, &c_b10, &t[i__ + mp * t_dim1], ldt, (ftnlen)1);

/*        B1 */

#line 290 "dtplqt2.f"
	i__2 = i__ - 1;
#line 290 "dtplqt2.f"
	i__3 = *n - *l;
#line 290 "dtplqt2.f"
	dgemv_("N", &i__2, &i__3, &alpha, &b[b_offset], ldb, &b[i__ + b_dim1],
		 ldb, &c_b4, &t[i__ + t_dim1], ldt, (ftnlen)1);

/*        T(1:I-1,I) := T(1:I-1,1:I-1) * T(I,1:I-1) */

#line 295 "dtplqt2.f"
	i__2 = i__ - 1;
#line 295 "dtplqt2.f"
	dtrmv_("L", "T", "N", &i__2, &t[t_offset], ldt, &t[i__ + t_dim1], ldt,
		 (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*        T(I,I) = tau(I) */

#line 299 "dtplqt2.f"
	t[i__ + i__ * t_dim1] = t[i__ * t_dim1 + 1];
#line 300 "dtplqt2.f"
	t[i__ * t_dim1 + 1] = 0.;
#line 301 "dtplqt2.f"
    }
#line 302 "dtplqt2.f"
    i__1 = *m;
#line 302 "dtplqt2.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 303 "dtplqt2.f"
	i__2 = *m;
#line 303 "dtplqt2.f"
	for (j = i__ + 1; j <= i__2; ++j) {
#line 304 "dtplqt2.f"
	    t[i__ + j * t_dim1] = t[j + i__ * t_dim1];
#line 305 "dtplqt2.f"
	    t[j + i__ * t_dim1] = 0.;
#line 306 "dtplqt2.f"
	}
#line 307 "dtplqt2.f"
    }

/*     End of DTPLQT2 */

#line 312 "dtplqt2.f"
    return 0;
} /* dtplqt2_ */

