#line 1 "cunmtr.f"
/* cunmtr.f -- translated by f2c (version 20100827).
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

#line 1 "cunmtr.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;

/* > \brief \b CUNMTR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CUNMTR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cunmtr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cunmtr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cunmtr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CUNMTR( SIDE, UPLO, TRANS, M, N, A, LDA, TAU, C, LDC, */
/*                          WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          SIDE, TRANS, UPLO */
/*       INTEGER            INFO, LDA, LDC, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            A( LDA, * ), C( LDC, * ), TAU( * ), */
/*      $                   WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CUNMTR overwrites the general complex M-by-N matrix C with */
/* > */
/* >                 SIDE = 'L'     SIDE = 'R' */
/* > TRANS = 'N':      Q * C          C * Q */
/* > TRANS = 'C':      Q**H * C       C * Q**H */
/* > */
/* > where Q is a complex unitary matrix of order nq, with nq = m if */
/* > SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of */
/* > nq-1 elementary reflectors, as returned by CHETRD: */
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
/* >                 from CHETRD; */
/* >          = 'L': Lower triangle of A contains elementary reflectors */
/* >                 from CHETRD. */
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
/* >          A is COMPLEX array, dimension */
/* >                               (LDA,M) if SIDE = 'L' */
/* >                               (LDA,N) if SIDE = 'R' */
/* >          The vectors which define the elementary reflectors, as */
/* >          returned by CHETRD. */
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
/* >          TAU is COMPLEX array, dimension */
/* >                               (M-1) if SIDE = 'L' */
/* >                               (N-1) if SIDE = 'R' */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i), as returned by CHETRD. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is COMPLEX array, dimension (LDC,N) */
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
/* >          WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
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

/* > \date December 2016 */

/* > \ingroup complexOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int cunmtr_(char *side, char *uplo, char *trans, integer *m, 
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
    extern /* Subroutine */ int cunmql_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen), cunmqr_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen);
    static integer lwkopt;
    static logical lquery;


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

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

#line 210 "cunmtr.f"
    /* Parameter adjustments */
#line 210 "cunmtr.f"
    a_dim1 = *lda;
#line 210 "cunmtr.f"
    a_offset = 1 + a_dim1;
#line 210 "cunmtr.f"
    a -= a_offset;
#line 210 "cunmtr.f"
    --tau;
#line 210 "cunmtr.f"
    c_dim1 = *ldc;
#line 210 "cunmtr.f"
    c_offset = 1 + c_dim1;
#line 210 "cunmtr.f"
    c__ -= c_offset;
#line 210 "cunmtr.f"
    --work;
#line 210 "cunmtr.f"

#line 210 "cunmtr.f"
    /* Function Body */
#line 210 "cunmtr.f"
    *info = 0;
#line 211 "cunmtr.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 212 "cunmtr.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 213 "cunmtr.f"
    lquery = *lwork == -1;

/*     NQ is the order of Q and NW is the minimum dimension of WORK */

#line 217 "cunmtr.f"
    if (left) {
#line 218 "cunmtr.f"
	nq = *m;
#line 219 "cunmtr.f"
	nw = *n;
#line 220 "cunmtr.f"
    } else {
#line 221 "cunmtr.f"
	nq = *n;
#line 222 "cunmtr.f"
	nw = *m;
#line 223 "cunmtr.f"
    }
#line 224 "cunmtr.f"
    if (! left && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 225 "cunmtr.f"
	*info = -1;
#line 226 "cunmtr.f"
    } else if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 227 "cunmtr.f"
	*info = -2;
#line 228 "cunmtr.f"
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "C", (ftnlen)1, (ftnlen)1)) {
#line 230 "cunmtr.f"
	*info = -3;
#line 231 "cunmtr.f"
    } else if (*m < 0) {
#line 232 "cunmtr.f"
	*info = -4;
#line 233 "cunmtr.f"
    } else if (*n < 0) {
#line 234 "cunmtr.f"
	*info = -5;
#line 235 "cunmtr.f"
    } else if (*lda < max(1,nq)) {
#line 236 "cunmtr.f"
	*info = -7;
#line 237 "cunmtr.f"
    } else if (*ldc < max(1,*m)) {
#line 238 "cunmtr.f"
	*info = -10;
#line 239 "cunmtr.f"
    } else if (*lwork < max(1,nw) && ! lquery) {
#line 240 "cunmtr.f"
	*info = -12;
#line 241 "cunmtr.f"
    }

#line 243 "cunmtr.f"
    if (*info == 0) {
#line 244 "cunmtr.f"
	if (upper) {
#line 245 "cunmtr.f"
	    if (left) {
/* Writing concatenation */
#line 246 "cunmtr.f"
		i__1[0] = 1, a__1[0] = side;
#line 246 "cunmtr.f"
		i__1[1] = 1, a__1[1] = trans;
#line 246 "cunmtr.f"
		s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
#line 246 "cunmtr.f"
		i__2 = *m - 1;
#line 246 "cunmtr.f"
		i__3 = *m - 1;
#line 246 "cunmtr.f"
		nb = ilaenv_(&c__1, "CUNMQL", ch__1, &i__2, n, &i__3, &c_n1, (
			ftnlen)6, (ftnlen)2);
#line 248 "cunmtr.f"
	    } else {
/* Writing concatenation */
#line 249 "cunmtr.f"
		i__1[0] = 1, a__1[0] = side;
#line 249 "cunmtr.f"
		i__1[1] = 1, a__1[1] = trans;
#line 249 "cunmtr.f"
		s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
#line 249 "cunmtr.f"
		i__2 = *n - 1;
#line 249 "cunmtr.f"
		i__3 = *n - 1;
#line 249 "cunmtr.f"
		nb = ilaenv_(&c__1, "CUNMQL", ch__1, m, &i__2, &i__3, &c_n1, (
			ftnlen)6, (ftnlen)2);
#line 251 "cunmtr.f"
	    }
#line 252 "cunmtr.f"
	} else {
#line 253 "cunmtr.f"
	    if (left) {
/* Writing concatenation */
#line 254 "cunmtr.f"
		i__1[0] = 1, a__1[0] = side;
#line 254 "cunmtr.f"
		i__1[1] = 1, a__1[1] = trans;
#line 254 "cunmtr.f"
		s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
#line 254 "cunmtr.f"
		i__2 = *m - 1;
#line 254 "cunmtr.f"
		i__3 = *m - 1;
#line 254 "cunmtr.f"
		nb = ilaenv_(&c__1, "CUNMQR", ch__1, &i__2, n, &i__3, &c_n1, (
			ftnlen)6, (ftnlen)2);
#line 256 "cunmtr.f"
	    } else {
/* Writing concatenation */
#line 257 "cunmtr.f"
		i__1[0] = 1, a__1[0] = side;
#line 257 "cunmtr.f"
		i__1[1] = 1, a__1[1] = trans;
#line 257 "cunmtr.f"
		s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
#line 257 "cunmtr.f"
		i__2 = *n - 1;
#line 257 "cunmtr.f"
		i__3 = *n - 1;
#line 257 "cunmtr.f"
		nb = ilaenv_(&c__1, "CUNMQR", ch__1, m, &i__2, &i__3, &c_n1, (
			ftnlen)6, (ftnlen)2);
#line 259 "cunmtr.f"
	    }
#line 260 "cunmtr.f"
	}
#line 261 "cunmtr.f"
	lwkopt = max(1,nw) * nb;
#line 262 "cunmtr.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 263 "cunmtr.f"
    }

#line 265 "cunmtr.f"
    if (*info != 0) {
#line 266 "cunmtr.f"
	i__2 = -(*info);
#line 266 "cunmtr.f"
	xerbla_("CUNMTR", &i__2, (ftnlen)6);
#line 267 "cunmtr.f"
	return 0;
#line 268 "cunmtr.f"
    } else if (lquery) {
#line 269 "cunmtr.f"
	return 0;
#line 270 "cunmtr.f"
    }

/*     Quick return if possible */

#line 274 "cunmtr.f"
    if (*m == 0 || *n == 0 || nq == 1) {
#line 275 "cunmtr.f"
	work[1].r = 1., work[1].i = 0.;
#line 276 "cunmtr.f"
	return 0;
#line 277 "cunmtr.f"
    }

#line 279 "cunmtr.f"
    if (left) {
#line 280 "cunmtr.f"
	mi = *m - 1;
#line 281 "cunmtr.f"
	ni = *n;
#line 282 "cunmtr.f"
    } else {
#line 283 "cunmtr.f"
	mi = *m;
#line 284 "cunmtr.f"
	ni = *n - 1;
#line 285 "cunmtr.f"
    }

#line 287 "cunmtr.f"
    if (upper) {

/*        Q was determined by a call to CHETRD with UPLO = 'U' */

#line 291 "cunmtr.f"
	i__2 = nq - 1;
#line 291 "cunmtr.f"
	cunmql_(side, trans, &mi, &ni, &i__2, &a[(a_dim1 << 1) + 1], lda, &
		tau[1], &c__[c_offset], ldc, &work[1], lwork, &iinfo, (ftnlen)
		1, (ftnlen)1);
#line 293 "cunmtr.f"
    } else {

/*        Q was determined by a call to CHETRD with UPLO = 'L' */

#line 297 "cunmtr.f"
	if (left) {
#line 298 "cunmtr.f"
	    i1 = 2;
#line 299 "cunmtr.f"
	    i2 = 1;
#line 300 "cunmtr.f"
	} else {
#line 301 "cunmtr.f"
	    i1 = 1;
#line 302 "cunmtr.f"
	    i2 = 2;
#line 303 "cunmtr.f"
	}
#line 304 "cunmtr.f"
	i__2 = nq - 1;
#line 304 "cunmtr.f"
	cunmqr_(side, trans, &mi, &ni, &i__2, &a[a_dim1 + 2], lda, &tau[1], &
		c__[i1 + i2 * c_dim1], ldc, &work[1], lwork, &iinfo, (ftnlen)
		1, (ftnlen)1);
#line 306 "cunmtr.f"
    }
#line 307 "cunmtr.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 308 "cunmtr.f"
    return 0;

/*     End of CUNMTR */

} /* cunmtr_ */

