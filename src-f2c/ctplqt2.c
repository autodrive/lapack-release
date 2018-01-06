#line 1 "ctplqt2.f"
/* ctplqt2.f -- translated by f2c (version 20100827).
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

#line 1 "ctplqt2.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CTPLQT2( M, N, L, A, LDA, B, LDB, T, LDT, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER   INFO, LDA, LDB, LDT, N, M, L */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX   A( LDA, * ), B( LDB, * ), T( LDT, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CTPLQT2 computes a LQ a factorization of a complex "triangular-pentagonal" */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
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
/* >          B is COMPLEX array, dimension (LDB,N) */
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
/* >          T is COMPLEX array, dimension (LDT,M) */
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
/* Subroutine */ int ctplqt2_(integer *m, integer *n, integer *l, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *t, integer *ldt, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, t_dim1, t_offset, i__1, i__2, 
	    i__3, i__4, i__5;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, p, mp, np;
    extern /* Subroutine */ int cgerc_(integer *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static doublecomplex alpha;
    extern /* Subroutine */ int cgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen), 
	    ctrmv_(char *, char *, char *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, ftnlen, ftnlen, ftnlen), 
	    clarfg_(integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *), xerbla_(char *, integer *, ftnlen);


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

#line 195 "ctplqt2.f"
    /* Parameter adjustments */
#line 195 "ctplqt2.f"
    a_dim1 = *lda;
#line 195 "ctplqt2.f"
    a_offset = 1 + a_dim1;
#line 195 "ctplqt2.f"
    a -= a_offset;
#line 195 "ctplqt2.f"
    b_dim1 = *ldb;
#line 195 "ctplqt2.f"
    b_offset = 1 + b_dim1;
#line 195 "ctplqt2.f"
    b -= b_offset;
#line 195 "ctplqt2.f"
    t_dim1 = *ldt;
#line 195 "ctplqt2.f"
    t_offset = 1 + t_dim1;
#line 195 "ctplqt2.f"
    t -= t_offset;
#line 195 "ctplqt2.f"

#line 195 "ctplqt2.f"
    /* Function Body */
#line 195 "ctplqt2.f"
    *info = 0;
#line 196 "ctplqt2.f"
    if (*m < 0) {
#line 197 "ctplqt2.f"
	*info = -1;
#line 198 "ctplqt2.f"
    } else if (*n < 0) {
#line 199 "ctplqt2.f"
	*info = -2;
#line 200 "ctplqt2.f"
    } else if (*l < 0 || *l > min(*m,*n)) {
#line 201 "ctplqt2.f"
	*info = -3;
#line 202 "ctplqt2.f"
    } else if (*lda < max(1,*m)) {
#line 203 "ctplqt2.f"
	*info = -5;
#line 204 "ctplqt2.f"
    } else if (*ldb < max(1,*m)) {
#line 205 "ctplqt2.f"
	*info = -7;
#line 206 "ctplqt2.f"
    } else if (*ldt < max(1,*m)) {
#line 207 "ctplqt2.f"
	*info = -9;
#line 208 "ctplqt2.f"
    }
#line 209 "ctplqt2.f"
    if (*info != 0) {
#line 210 "ctplqt2.f"
	i__1 = -(*info);
#line 210 "ctplqt2.f"
	xerbla_("CTPLQT2", &i__1, (ftnlen)7);
#line 211 "ctplqt2.f"
	return 0;
#line 212 "ctplqt2.f"
    }

/*     Quick return if possible */

#line 216 "ctplqt2.f"
    if (*n == 0 || *m == 0) {
#line 216 "ctplqt2.f"
	return 0;
#line 216 "ctplqt2.f"
    }

#line 218 "ctplqt2.f"
    i__1 = *m;
#line 218 "ctplqt2.f"
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        Generate elementary reflector H(I) to annihilate B(I,:) */

#line 222 "ctplqt2.f"
	p = *n - *l + min(*l,i__);
#line 223 "ctplqt2.f"
	i__2 = p + 1;
#line 223 "ctplqt2.f"
	clarfg_(&i__2, &a[i__ + i__ * a_dim1], &b[i__ + b_dim1], ldb, &t[i__ *
		 t_dim1 + 1]);
#line 224 "ctplqt2.f"
	i__2 = i__ * t_dim1 + 1;
#line 224 "ctplqt2.f"
	d_cnjg(&z__1, &t[i__ * t_dim1 + 1]);
#line 224 "ctplqt2.f"
	t[i__2].r = z__1.r, t[i__2].i = z__1.i;
#line 225 "ctplqt2.f"
	if (i__ < *m) {
#line 226 "ctplqt2.f"
	    i__2 = p;
#line 226 "ctplqt2.f"
	    for (j = 1; j <= i__2; ++j) {
#line 227 "ctplqt2.f"
		i__3 = i__ + j * b_dim1;
#line 227 "ctplqt2.f"
		d_cnjg(&z__1, &b[i__ + j * b_dim1]);
#line 227 "ctplqt2.f"
		b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 228 "ctplqt2.f"
	    }

/*           W(M-I:1) := C(I+1:M,I:N) * C(I,I:N) [use W = T(M,:)] */

#line 232 "ctplqt2.f"
	    i__2 = *m - i__;
#line 232 "ctplqt2.f"
	    for (j = 1; j <= i__2; ++j) {
#line 233 "ctplqt2.f"
		i__3 = *m + j * t_dim1;
#line 233 "ctplqt2.f"
		i__4 = i__ + j + i__ * a_dim1;
#line 233 "ctplqt2.f"
		t[i__3].r = a[i__4].r, t[i__3].i = a[i__4].i;
#line 234 "ctplqt2.f"
	    }
#line 235 "ctplqt2.f"
	    i__2 = *m - i__;
#line 235 "ctplqt2.f"
	    cgemv_("N", &i__2, &p, &c_b2, &b[i__ + 1 + b_dim1], ldb, &b[i__ + 
		    b_dim1], ldb, &c_b2, &t[*m + t_dim1], ldt, (ftnlen)1);

/*           C(I+1:M,I:N) = C(I+1:M,I:N) + alpha * C(I,I:N)*W(M-1:1)^H */

#line 240 "ctplqt2.f"
	    i__2 = i__ * t_dim1 + 1;
#line 240 "ctplqt2.f"
	    z__1.r = -t[i__2].r, z__1.i = -t[i__2].i;
#line 240 "ctplqt2.f"
	    alpha.r = z__1.r, alpha.i = z__1.i;
#line 241 "ctplqt2.f"
	    i__2 = *m - i__;
#line 241 "ctplqt2.f"
	    for (j = 1; j <= i__2; ++j) {
#line 242 "ctplqt2.f"
		i__3 = i__ + j + i__ * a_dim1;
#line 242 "ctplqt2.f"
		i__4 = i__ + j + i__ * a_dim1;
#line 242 "ctplqt2.f"
		i__5 = *m + j * t_dim1;
#line 242 "ctplqt2.f"
		z__2.r = alpha.r * t[i__5].r - alpha.i * t[i__5].i, z__2.i = 
			alpha.r * t[i__5].i + alpha.i * t[i__5].r;
#line 242 "ctplqt2.f"
		z__1.r = a[i__4].r + z__2.r, z__1.i = a[i__4].i + z__2.i;
#line 242 "ctplqt2.f"
		a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 243 "ctplqt2.f"
	    }
#line 244 "ctplqt2.f"
	    i__2 = *m - i__;
#line 244 "ctplqt2.f"
	    z__1.r = alpha.r, z__1.i = alpha.i;
#line 244 "ctplqt2.f"
	    cgerc_(&i__2, &p, &z__1, &t[*m + t_dim1], ldt, &b[i__ + b_dim1], 
		    ldb, &b[i__ + 1 + b_dim1], ldb);
#line 246 "ctplqt2.f"
	    i__2 = p;
#line 246 "ctplqt2.f"
	    for (j = 1; j <= i__2; ++j) {
#line 247 "ctplqt2.f"
		i__3 = i__ + j * b_dim1;
#line 247 "ctplqt2.f"
		d_cnjg(&z__1, &b[i__ + j * b_dim1]);
#line 247 "ctplqt2.f"
		b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 248 "ctplqt2.f"
	    }
#line 249 "ctplqt2.f"
	}
#line 250 "ctplqt2.f"
    }

#line 252 "ctplqt2.f"
    i__1 = *m;
#line 252 "ctplqt2.f"
    for (i__ = 2; i__ <= i__1; ++i__) {

/*        T(I,1:I-1) := C(I:I-1,1:N)**H * (alpha * C(I,I:N)) */

#line 256 "ctplqt2.f"
	i__2 = i__ * t_dim1 + 1;
#line 256 "ctplqt2.f"
	z__1.r = -t[i__2].r, z__1.i = -t[i__2].i;
#line 256 "ctplqt2.f"
	alpha.r = z__1.r, alpha.i = z__1.i;
#line 257 "ctplqt2.f"
	i__2 = i__ - 1;
#line 257 "ctplqt2.f"
	for (j = 1; j <= i__2; ++j) {
#line 258 "ctplqt2.f"
	    i__3 = i__ + j * t_dim1;
#line 258 "ctplqt2.f"
	    t[i__3].r = 0., t[i__3].i = 0.;
#line 259 "ctplqt2.f"
	}
/* Computing MIN */
#line 260 "ctplqt2.f"
	i__2 = i__ - 1;
#line 260 "ctplqt2.f"
	p = min(i__2,*l);
/* Computing MIN */
#line 261 "ctplqt2.f"
	i__2 = *n - *l + 1;
#line 261 "ctplqt2.f"
	np = min(i__2,*n);
/* Computing MIN */
#line 262 "ctplqt2.f"
	i__2 = p + 1;
#line 262 "ctplqt2.f"
	mp = min(i__2,*m);
#line 263 "ctplqt2.f"
	i__2 = *n - *l + p;
#line 263 "ctplqt2.f"
	for (j = 1; j <= i__2; ++j) {
#line 264 "ctplqt2.f"
	    i__3 = i__ + j * b_dim1;
#line 264 "ctplqt2.f"
	    d_cnjg(&z__1, &b[i__ + j * b_dim1]);
#line 264 "ctplqt2.f"
	    b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 265 "ctplqt2.f"
	}

/*        Triangular part of B2 */

#line 269 "ctplqt2.f"
	i__2 = p;
#line 269 "ctplqt2.f"
	for (j = 1; j <= i__2; ++j) {
#line 270 "ctplqt2.f"
	    i__3 = i__ + j * t_dim1;
#line 270 "ctplqt2.f"
	    i__4 = i__ + (*n - *l + j) * b_dim1;
#line 270 "ctplqt2.f"
	    z__1.r = alpha.r * b[i__4].r - alpha.i * b[i__4].i, z__1.i = 
		    alpha.r * b[i__4].i + alpha.i * b[i__4].r;
#line 270 "ctplqt2.f"
	    t[i__3].r = z__1.r, t[i__3].i = z__1.i;
#line 271 "ctplqt2.f"
	}
#line 272 "ctplqt2.f"
	ctrmv_("L", "N", "N", &p, &b[np * b_dim1 + 1], ldb, &t[i__ + t_dim1], 
		ldt, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*        Rectangular part of B2 */

#line 277 "ctplqt2.f"
	i__2 = i__ - 1 - p;
#line 277 "ctplqt2.f"
	cgemv_("N", &i__2, l, &alpha, &b[mp + np * b_dim1], ldb, &b[i__ + np *
		 b_dim1], ldb, &c_b1, &t[i__ + mp * t_dim1], ldt, (ftnlen)1);

/*        B1 */

#line 283 "ctplqt2.f"
	i__2 = i__ - 1;
#line 283 "ctplqt2.f"
	i__3 = *n - *l;
#line 283 "ctplqt2.f"
	cgemv_("N", &i__2, &i__3, &alpha, &b[b_offset], ldb, &b[i__ + b_dim1],
		 ldb, &c_b2, &t[i__ + t_dim1], ldt, (ftnlen)1);


/*        T(1:I-1,I) := T(1:I-1,1:I-1) * T(I,1:I-1) */

#line 290 "ctplqt2.f"
	i__2 = i__ - 1;
#line 290 "ctplqt2.f"
	for (j = 1; j <= i__2; ++j) {
#line 291 "ctplqt2.f"
	    i__3 = i__ + j * t_dim1;
#line 291 "ctplqt2.f"
	    d_cnjg(&z__1, &t[i__ + j * t_dim1]);
#line 291 "ctplqt2.f"
	    t[i__3].r = z__1.r, t[i__3].i = z__1.i;
#line 292 "ctplqt2.f"
	}
#line 293 "ctplqt2.f"
	i__2 = i__ - 1;
#line 293 "ctplqt2.f"
	ctrmv_("L", "C", "N", &i__2, &t[t_offset], ldt, &t[i__ + t_dim1], ldt,
		 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 294 "ctplqt2.f"
	i__2 = i__ - 1;
#line 294 "ctplqt2.f"
	for (j = 1; j <= i__2; ++j) {
#line 295 "ctplqt2.f"
	    i__3 = i__ + j * t_dim1;
#line 295 "ctplqt2.f"
	    d_cnjg(&z__1, &t[i__ + j * t_dim1]);
#line 295 "ctplqt2.f"
	    t[i__3].r = z__1.r, t[i__3].i = z__1.i;
#line 296 "ctplqt2.f"
	}
#line 297 "ctplqt2.f"
	i__2 = *n - *l + p;
#line 297 "ctplqt2.f"
	for (j = 1; j <= i__2; ++j) {
#line 298 "ctplqt2.f"
	    i__3 = i__ + j * b_dim1;
#line 298 "ctplqt2.f"
	    d_cnjg(&z__1, &b[i__ + j * b_dim1]);
#line 298 "ctplqt2.f"
	    b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 299 "ctplqt2.f"
	}

/*        T(I,I) = tau(I) */

#line 303 "ctplqt2.f"
	i__2 = i__ + i__ * t_dim1;
#line 303 "ctplqt2.f"
	i__3 = i__ * t_dim1 + 1;
#line 303 "ctplqt2.f"
	t[i__2].r = t[i__3].r, t[i__2].i = t[i__3].i;
#line 304 "ctplqt2.f"
	i__2 = i__ * t_dim1 + 1;
#line 304 "ctplqt2.f"
	t[i__2].r = 0., t[i__2].i = 0.;
#line 305 "ctplqt2.f"
    }
#line 306 "ctplqt2.f"
    i__1 = *m;
#line 306 "ctplqt2.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 307 "ctplqt2.f"
	i__2 = *m;
#line 307 "ctplqt2.f"
	for (j = i__ + 1; j <= i__2; ++j) {
#line 308 "ctplqt2.f"
	    i__3 = i__ + j * t_dim1;
#line 308 "ctplqt2.f"
	    i__4 = j + i__ * t_dim1;
#line 308 "ctplqt2.f"
	    t[i__3].r = t[i__4].r, t[i__3].i = t[i__4].i;
#line 309 "ctplqt2.f"
	    i__3 = j + i__ * t_dim1;
#line 309 "ctplqt2.f"
	    t[i__3].r = 0., t[i__3].i = 0.;
#line 310 "ctplqt2.f"
	}
#line 311 "ctplqt2.f"
    }

/*     End of CTPLQT2 */

#line 316 "ctplqt2.f"
    return 0;
} /* ctplqt2_ */

