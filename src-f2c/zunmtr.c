#line 1 "zunmtr.f"
/* zunmtr.f -- translated by f2c (version 20100827).
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

#line 1 "zunmtr.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;

/* > \brief \b ZUNMTR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZUNMTR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zunmtr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zunmtr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zunmtr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZUNMTR( SIDE, UPLO, TRANS, M, N, A, LDA, TAU, C, LDC, */
/*                          WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          SIDE, TRANS, UPLO */
/*       INTEGER            INFO, LDA, LDC, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZUNMTR overwrites the general complex M-by-N matrix C with */
/* > */
/* >                 SIDE = 'L'     SIDE = 'R' */
/* > TRANS = 'N':      Q * C          C * Q */
/* > TRANS = 'C':      Q**H * C       C * Q**H */
/* > */
/* > where Q is a complex unitary matrix of order nq, with nq = m if */
/* > SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of */
/* > nq-1 elementary reflectors, as returned by ZHETRD: */
/* > */
/* > if UPLO = 'U', Q = H(nq-1) . . . H(2) H(1); */
/* > */
/* > if UPLO = 'L', Q = H(1) H(2) . . . H(nq-1). */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >          = 'L': apply Q or Q**H from the Left; */
/* >          = 'R': apply Q or Q**H from the Right. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U': Upper triangle of A contains elementary reflectors */
/* >                 from ZHETRD; */
/* >          = 'L': Lower triangle of A contains elementary reflectors */
/* >                 from ZHETRD. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          = 'N':  No transpose, apply Q; */
/* >          = 'C':  Conjugate transpose, apply Q**H. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix C. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix C. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension */
/* >                               (LDA,M) if SIDE = 'L' */
/* >                               (LDA,N) if SIDE = 'R' */
/* >          The vectors which define the elementary reflectors, as */
/* >          returned by ZHETRD. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A. */
/* >          LDA >= max(1,M) if SIDE = 'L'; LDA >= max(1,N) if SIDE = 'R'. */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* >          TAU is COMPLEX*16 array, dimension */
/* >                               (M-1) if SIDE = 'L' */
/* >                               (N-1) if SIDE = 'R' */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i), as returned by ZHETRD. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is COMPLEX*16 array, dimension (LDC,N) */
/* >          On entry, the M-by-N matrix C. */
/* >          On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* >          LDC is INTEGER */
/* >          The leading dimension of the array C. LDC >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. */
/* >          If SIDE = 'L', LWORK >= max(1,N); */
/* >          if SIDE = 'R', LWORK >= max(1,M). */
/* >          For optimum performance LWORK >= N*NB if SIDE = 'L', and */
/* >          LWORK >=M*NB if SIDE = 'R', where NB is the optimal */
/* >          blocksize. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup complex16OTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int zunmtr_(char *side, char *uplo, char *trans, integer *m, 
	integer *n, doublecomplex *a, integer *lda, doublecomplex *tau, 
	doublecomplex *c__, integer *ldc, doublecomplex *work, integer *lwork,
	 integer *info, ftnlen side_len, ftnlen uplo_len, ftnlen trans_len)
{
    /* System generated locals */
    address a__1[2];
    integer a_dim1, a_offset, c_dim1, c_offset, i__1[2], i__2, i__3;
    char ch__1[2];

    /* Builtin functions */
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);

    /* Local variables */
    static integer i1, i2, nb, mi, ni, nq, nw;
    static logical left;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer lwkopt;
    static logical lquery;
    extern /* Subroutine */ int zunmql_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen), zunmqr_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen);


/*  -- LAPACK computational routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments */

#line 208 "zunmtr.f"
    /* Parameter adjustments */
#line 208 "zunmtr.f"
    a_dim1 = *lda;
#line 208 "zunmtr.f"
    a_offset = 1 + a_dim1;
#line 208 "zunmtr.f"
    a -= a_offset;
#line 208 "zunmtr.f"
    --tau;
#line 208 "zunmtr.f"
    c_dim1 = *ldc;
#line 208 "zunmtr.f"
    c_offset = 1 + c_dim1;
#line 208 "zunmtr.f"
    c__ -= c_offset;
#line 208 "zunmtr.f"
    --work;
#line 208 "zunmtr.f"

#line 208 "zunmtr.f"
    /* Function Body */
#line 208 "zunmtr.f"
    *info = 0;
#line 209 "zunmtr.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 210 "zunmtr.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 211 "zunmtr.f"
    lquery = *lwork == -1;

/*     NQ is the order of Q and NW is the minimum dimension of WORK */

#line 215 "zunmtr.f"
    if (left) {
#line 216 "zunmtr.f"
	nq = *m;
#line 217 "zunmtr.f"
	nw = *n;
#line 218 "zunmtr.f"
    } else {
#line 219 "zunmtr.f"
	nq = *n;
#line 220 "zunmtr.f"
	nw = *m;
#line 221 "zunmtr.f"
    }
#line 222 "zunmtr.f"
    if (! left && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 223 "zunmtr.f"
	*info = -1;
#line 224 "zunmtr.f"
    } else if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 225 "zunmtr.f"
	*info = -2;
#line 226 "zunmtr.f"
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "C", (ftnlen)1, (ftnlen)1)) {
#line 228 "zunmtr.f"
	*info = -3;
#line 229 "zunmtr.f"
    } else if (*m < 0) {
#line 230 "zunmtr.f"
	*info = -4;
#line 231 "zunmtr.f"
    } else if (*n < 0) {
#line 232 "zunmtr.f"
	*info = -5;
#line 233 "zunmtr.f"
    } else if (*lda < max(1,nq)) {
#line 234 "zunmtr.f"
	*info = -7;
#line 235 "zunmtr.f"
    } else if (*ldc < max(1,*m)) {
#line 236 "zunmtr.f"
	*info = -10;
#line 237 "zunmtr.f"
    } else if (*lwork < max(1,nw) && ! lquery) {
#line 238 "zunmtr.f"
	*info = -12;
#line 239 "zunmtr.f"
    }

#line 241 "zunmtr.f"
    if (*info == 0) {
#line 242 "zunmtr.f"
	if (upper) {
#line 243 "zunmtr.f"
	    if (left) {
/* Writing concatenation */
#line 244 "zunmtr.f"
		i__1[0] = 1, a__1[0] = side;
#line 244 "zunmtr.f"
		i__1[1] = 1, a__1[1] = trans;
#line 244 "zunmtr.f"
		s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
#line 244 "zunmtr.f"
		i__2 = *m - 1;
#line 244 "zunmtr.f"
		i__3 = *m - 1;
#line 244 "zunmtr.f"
		nb = ilaenv_(&c__1, "ZUNMQL", ch__1, &i__2, n, &i__3, &c_n1, (
			ftnlen)6, (ftnlen)2);
#line 246 "zunmtr.f"
	    } else {
/* Writing concatenation */
#line 247 "zunmtr.f"
		i__1[0] = 1, a__1[0] = side;
#line 247 "zunmtr.f"
		i__1[1] = 1, a__1[1] = trans;
#line 247 "zunmtr.f"
		s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
#line 247 "zunmtr.f"
		i__2 = *n - 1;
#line 247 "zunmtr.f"
		i__3 = *n - 1;
#line 247 "zunmtr.f"
		nb = ilaenv_(&c__1, "ZUNMQL", ch__1, m, &i__2, &i__3, &c_n1, (
			ftnlen)6, (ftnlen)2);
#line 249 "zunmtr.f"
	    }
#line 250 "zunmtr.f"
	} else {
#line 251 "zunmtr.f"
	    if (left) {
/* Writing concatenation */
#line 252 "zunmtr.f"
		i__1[0] = 1, a__1[0] = side;
#line 252 "zunmtr.f"
		i__1[1] = 1, a__1[1] = trans;
#line 252 "zunmtr.f"
		s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
#line 252 "zunmtr.f"
		i__2 = *m - 1;
#line 252 "zunmtr.f"
		i__3 = *m - 1;
#line 252 "zunmtr.f"
		nb = ilaenv_(&c__1, "ZUNMQR", ch__1, &i__2, n, &i__3, &c_n1, (
			ftnlen)6, (ftnlen)2);
#line 254 "zunmtr.f"
	    } else {
/* Writing concatenation */
#line 255 "zunmtr.f"
		i__1[0] = 1, a__1[0] = side;
#line 255 "zunmtr.f"
		i__1[1] = 1, a__1[1] = trans;
#line 255 "zunmtr.f"
		s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
#line 255 "zunmtr.f"
		i__2 = *n - 1;
#line 255 "zunmtr.f"
		i__3 = *n - 1;
#line 255 "zunmtr.f"
		nb = ilaenv_(&c__1, "ZUNMQR", ch__1, m, &i__2, &i__3, &c_n1, (
			ftnlen)6, (ftnlen)2);
#line 257 "zunmtr.f"
	    }
#line 258 "zunmtr.f"
	}
#line 259 "zunmtr.f"
	lwkopt = max(1,nw) * nb;
#line 260 "zunmtr.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 261 "zunmtr.f"
    }

#line 263 "zunmtr.f"
    if (*info != 0) {
#line 264 "zunmtr.f"
	i__2 = -(*info);
#line 264 "zunmtr.f"
	xerbla_("ZUNMTR", &i__2, (ftnlen)6);
#line 265 "zunmtr.f"
	return 0;
#line 266 "zunmtr.f"
    } else if (lquery) {
#line 267 "zunmtr.f"
	return 0;
#line 268 "zunmtr.f"
    }

/*     Quick return if possible */

#line 272 "zunmtr.f"
    if (*m == 0 || *n == 0 || nq == 1) {
#line 273 "zunmtr.f"
	work[1].r = 1., work[1].i = 0.;
#line 274 "zunmtr.f"
	return 0;
#line 275 "zunmtr.f"
    }

#line 277 "zunmtr.f"
    if (left) {
#line 278 "zunmtr.f"
	mi = *m - 1;
#line 279 "zunmtr.f"
	ni = *n;
#line 280 "zunmtr.f"
    } else {
#line 281 "zunmtr.f"
	mi = *m;
#line 282 "zunmtr.f"
	ni = *n - 1;
#line 283 "zunmtr.f"
    }

#line 285 "zunmtr.f"
    if (upper) {

/*        Q was determined by a call to ZHETRD with UPLO = 'U' */

#line 289 "zunmtr.f"
	i__2 = nq - 1;
#line 289 "zunmtr.f"
	zunmql_(side, trans, &mi, &ni, &i__2, &a[(a_dim1 << 1) + 1], lda, &
		tau[1], &c__[c_offset], ldc, &work[1], lwork, &iinfo, (ftnlen)
		1, (ftnlen)1);
#line 291 "zunmtr.f"
    } else {

/*        Q was determined by a call to ZHETRD with UPLO = 'L' */

#line 295 "zunmtr.f"
	if (left) {
#line 296 "zunmtr.f"
	    i1 = 2;
#line 297 "zunmtr.f"
	    i2 = 1;
#line 298 "zunmtr.f"
	} else {
#line 299 "zunmtr.f"
	    i1 = 1;
#line 300 "zunmtr.f"
	    i2 = 2;
#line 301 "zunmtr.f"
	}
#line 302 "zunmtr.f"
	i__2 = nq - 1;
#line 302 "zunmtr.f"
	zunmqr_(side, trans, &mi, &ni, &i__2, &a[a_dim1 + 2], lda, &tau[1], &
		c__[i1 + i2 * c_dim1], ldc, &work[1], lwork, &iinfo, (ftnlen)
		1, (ftnlen)1);
#line 304 "zunmtr.f"
    }
#line 305 "zunmtr.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 306 "zunmtr.f"
    return 0;

/*     End of ZUNMTR */

} /* zunmtr_ */

