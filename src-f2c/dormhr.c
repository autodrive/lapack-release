#line 1 "dormhr.f"
/* dormhr.f -- translated by f2c (version 20100827).
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

#line 1 "dormhr.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;

/* > \brief \b DORMHR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DORMHR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dormhr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dormhr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dormhr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DORMHR( SIDE, TRANS, M, N, ILO, IHI, A, LDA, TAU, C, */
/*                          LDC, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          SIDE, TRANS */
/*       INTEGER            IHI, ILO, INFO, LDA, LDC, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DORMHR overwrites the general real M-by-N matrix C with */
/* > */
/* >                 SIDE = 'L'     SIDE = 'R' */
/* > TRANS = 'N':      Q * C          C * Q */
/* > TRANS = 'T':      Q**T * C       C * Q**T */
/* > */
/* > where Q is a real orthogonal matrix of order nq, with nq = m if */
/* > SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of */
/* > IHI-ILO elementary reflectors, as returned by DGEHRD: */
/* > */
/* > Q = H(ilo) H(ilo+1) . . . H(ihi-1). */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >          = 'L': apply Q or Q**T from the Left; */
/* >          = 'R': apply Q or Q**T from the Right. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          = 'N':  No transpose, apply Q; */
/* >          = 'T':  Transpose, apply Q**T. */
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
/* > \param[in] ILO */
/* > \verbatim */
/* >          ILO is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] IHI */
/* > \verbatim */
/* >          IHI is INTEGER */
/* > */
/* >          ILO and IHI must have the same values as in the previous call */
/* >          of DGEHRD. Q is equal to the unit matrix except in the */
/* >          submatrix Q(ilo+1:ihi,ilo+1:ihi). */
/* >          If SIDE = 'L', then 1 <= ILO <= IHI <= M, if M > 0, and */
/* >          ILO = 1 and IHI = 0, if M = 0; */
/* >          if SIDE = 'R', then 1 <= ILO <= IHI <= N, if N > 0, and */
/* >          ILO = 1 and IHI = 0, if N = 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension */
/* >                               (LDA,M) if SIDE = 'L' */
/* >                               (LDA,N) if SIDE = 'R' */
/* >          The vectors which define the elementary reflectors, as */
/* >          returned by DGEHRD. */
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
/* >          TAU is DOUBLE PRECISION array, dimension */
/* >                               (M-1) if SIDE = 'L' */
/* >                               (N-1) if SIDE = 'R' */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i), as returned by DGEHRD. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION array, dimension (LDC,N) */
/* >          On entry, the M-by-N matrix C. */
/* >          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q. */
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
/* >          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
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
/* >          LWORK >= M*NB if SIDE = 'R', where NB is the optimal */
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

/* > \ingroup doubleOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int dormhr_(char *side, char *trans, integer *m, integer *n, 
	integer *ilo, integer *ihi, doublereal *a, integer *lda, doublereal *
	tau, doublereal *c__, integer *ldc, doublereal *work, integer *lwork, 
	integer *info, ftnlen side_len, ftnlen trans_len)
{
    /* System generated locals */
    address a__1[2];
    integer a_dim1, a_offset, c_dim1, c_offset, i__1[2], i__2;
    char ch__1[2];

    /* Builtin functions */
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);

    /* Local variables */
    static integer i1, i2, nb, mi, nh, ni, nq, nw;
    static logical left;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int dormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static integer lwkopt;
    static logical lquery;


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

#line 215 "dormhr.f"
    /* Parameter adjustments */
#line 215 "dormhr.f"
    a_dim1 = *lda;
#line 215 "dormhr.f"
    a_offset = 1 + a_dim1;
#line 215 "dormhr.f"
    a -= a_offset;
#line 215 "dormhr.f"
    --tau;
#line 215 "dormhr.f"
    c_dim1 = *ldc;
#line 215 "dormhr.f"
    c_offset = 1 + c_dim1;
#line 215 "dormhr.f"
    c__ -= c_offset;
#line 215 "dormhr.f"
    --work;
#line 215 "dormhr.f"

#line 215 "dormhr.f"
    /* Function Body */
#line 215 "dormhr.f"
    *info = 0;
#line 216 "dormhr.f"
    nh = *ihi - *ilo;
#line 217 "dormhr.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 218 "dormhr.f"
    lquery = *lwork == -1;

/*     NQ is the order of Q and NW is the minimum dimension of WORK */

#line 222 "dormhr.f"
    if (left) {
#line 223 "dormhr.f"
	nq = *m;
#line 224 "dormhr.f"
	nw = *n;
#line 225 "dormhr.f"
    } else {
#line 226 "dormhr.f"
	nq = *n;
#line 227 "dormhr.f"
	nw = *m;
#line 228 "dormhr.f"
    }
#line 229 "dormhr.f"
    if (! left && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 230 "dormhr.f"
	*info = -1;
#line 231 "dormhr.f"
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "T", (ftnlen)1, (ftnlen)1)) {
#line 233 "dormhr.f"
	*info = -2;
#line 234 "dormhr.f"
    } else if (*m < 0) {
#line 235 "dormhr.f"
	*info = -3;
#line 236 "dormhr.f"
    } else if (*n < 0) {
#line 237 "dormhr.f"
	*info = -4;
#line 238 "dormhr.f"
    } else if (*ilo < 1 || *ilo > max(1,nq)) {
#line 239 "dormhr.f"
	*info = -5;
#line 240 "dormhr.f"
    } else if (*ihi < min(*ilo,nq) || *ihi > nq) {
#line 241 "dormhr.f"
	*info = -6;
#line 242 "dormhr.f"
    } else if (*lda < max(1,nq)) {
#line 243 "dormhr.f"
	*info = -8;
#line 244 "dormhr.f"
    } else if (*ldc < max(1,*m)) {
#line 245 "dormhr.f"
	*info = -11;
#line 246 "dormhr.f"
    } else if (*lwork < max(1,nw) && ! lquery) {
#line 247 "dormhr.f"
	*info = -13;
#line 248 "dormhr.f"
    }

#line 250 "dormhr.f"
    if (*info == 0) {
#line 251 "dormhr.f"
	if (left) {
/* Writing concatenation */
#line 252 "dormhr.f"
	    i__1[0] = 1, a__1[0] = side;
#line 252 "dormhr.f"
	    i__1[1] = 1, a__1[1] = trans;
#line 252 "dormhr.f"
	    s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
#line 252 "dormhr.f"
	    nb = ilaenv_(&c__1, "DORMQR", ch__1, &nh, n, &nh, &c_n1, (ftnlen)
		    6, (ftnlen)2);
#line 253 "dormhr.f"
	} else {
/* Writing concatenation */
#line 254 "dormhr.f"
	    i__1[0] = 1, a__1[0] = side;
#line 254 "dormhr.f"
	    i__1[1] = 1, a__1[1] = trans;
#line 254 "dormhr.f"
	    s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
#line 254 "dormhr.f"
	    nb = ilaenv_(&c__1, "DORMQR", ch__1, m, &nh, &nh, &c_n1, (ftnlen)
		    6, (ftnlen)2);
#line 255 "dormhr.f"
	}
#line 256 "dormhr.f"
	lwkopt = max(1,nw) * nb;
#line 257 "dormhr.f"
	work[1] = (doublereal) lwkopt;
#line 258 "dormhr.f"
    }

#line 260 "dormhr.f"
    if (*info != 0) {
#line 261 "dormhr.f"
	i__2 = -(*info);
#line 261 "dormhr.f"
	xerbla_("DORMHR", &i__2, (ftnlen)6);
#line 262 "dormhr.f"
	return 0;
#line 263 "dormhr.f"
    } else if (lquery) {
#line 264 "dormhr.f"
	return 0;
#line 265 "dormhr.f"
    }

/*     Quick return if possible */

#line 269 "dormhr.f"
    if (*m == 0 || *n == 0 || nh == 0) {
#line 270 "dormhr.f"
	work[1] = 1.;
#line 271 "dormhr.f"
	return 0;
#line 272 "dormhr.f"
    }

#line 274 "dormhr.f"
    if (left) {
#line 275 "dormhr.f"
	mi = nh;
#line 276 "dormhr.f"
	ni = *n;
#line 277 "dormhr.f"
	i1 = *ilo + 1;
#line 278 "dormhr.f"
	i2 = 1;
#line 279 "dormhr.f"
    } else {
#line 280 "dormhr.f"
	mi = *m;
#line 281 "dormhr.f"
	ni = nh;
#line 282 "dormhr.f"
	i1 = 1;
#line 283 "dormhr.f"
	i2 = *ilo + 1;
#line 284 "dormhr.f"
    }

#line 286 "dormhr.f"
    dormqr_(side, trans, &mi, &ni, &nh, &a[*ilo + 1 + *ilo * a_dim1], lda, &
	    tau[*ilo], &c__[i1 + i2 * c_dim1], ldc, &work[1], lwork, &iinfo, (
	    ftnlen)1, (ftnlen)1);

#line 289 "dormhr.f"
    work[1] = (doublereal) lwkopt;
#line 290 "dormhr.f"
    return 0;

/*     End of DORMHR */

} /* dormhr_ */

