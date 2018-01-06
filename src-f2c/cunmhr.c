#line 1 "cunmhr.f"
/* cunmhr.f -- translated by f2c (version 20100827).
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

#line 1 "cunmhr.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;

/* > \brief \b CUNMHR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CUNMHR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cunmhr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cunmhr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cunmhr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CUNMHR( SIDE, TRANS, M, N, ILO, IHI, A, LDA, TAU, C, */
/*                          LDC, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          SIDE, TRANS */
/*       INTEGER            IHI, ILO, INFO, LDA, LDC, LWORK, M, N */
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
/* > CUNMHR overwrites the general complex M-by-N matrix C with */
/* > */
/* >                 SIDE = 'L'     SIDE = 'R' */
/* > TRANS = 'N':      Q * C          C * Q */
/* > TRANS = 'C':      Q**H * C       C * Q**H */
/* > */
/* > where Q is a complex unitary matrix of order nq, with nq = m if */
/* > SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of */
/* > IHI-ILO elementary reflectors, as returned by CGEHRD: */
/* > */
/* > Q = H(ilo) H(ilo+1) . . . H(ihi-1). */
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
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          = 'N': apply Q  (No transpose) */
/* >          = 'C': apply Q**H (Conjugate transpose) */
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
/* >          of CGEHRD. Q is equal to the unit matrix except in the */
/* >          submatrix Q(ilo+1:ihi,ilo+1:ihi). */
/* >          If SIDE = 'L', then 1 <= ILO <= IHI <= M, if M > 0, and */
/* >          ILO = 1 and IHI = 0, if M = 0; */
/* >          if SIDE = 'R', then 1 <= ILO <= IHI <= N, if N > 0, and */
/* >          ILO = 1 and IHI = 0, if N = 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension */
/* >                               (LDA,M) if SIDE = 'L' */
/* >                               (LDA,N) if SIDE = 'R' */
/* >          The vectors which define the elementary reflectors, as */
/* >          returned by CGEHRD. */
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
/* >          reflector H(i), as returned by CGEHRD. */
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

/* > \ingroup complexOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int cunmhr_(char *side, char *trans, integer *m, integer *n, 
	integer *ilo, integer *ihi, doublecomplex *a, integer *lda, 
	doublecomplex *tau, doublecomplex *c__, integer *ldc, doublecomplex *
	work, integer *lwork, integer *info, ftnlen side_len, ftnlen 
	trans_len)
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
    extern /* Subroutine */ int cunmqr_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen);
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

#line 217 "cunmhr.f"
    /* Parameter adjustments */
#line 217 "cunmhr.f"
    a_dim1 = *lda;
#line 217 "cunmhr.f"
    a_offset = 1 + a_dim1;
#line 217 "cunmhr.f"
    a -= a_offset;
#line 217 "cunmhr.f"
    --tau;
#line 217 "cunmhr.f"
    c_dim1 = *ldc;
#line 217 "cunmhr.f"
    c_offset = 1 + c_dim1;
#line 217 "cunmhr.f"
    c__ -= c_offset;
#line 217 "cunmhr.f"
    --work;
#line 217 "cunmhr.f"

#line 217 "cunmhr.f"
    /* Function Body */
#line 217 "cunmhr.f"
    *info = 0;
#line 218 "cunmhr.f"
    nh = *ihi - *ilo;
#line 219 "cunmhr.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 220 "cunmhr.f"
    lquery = *lwork == -1;

/*     NQ is the order of Q and NW is the minimum dimension of WORK */

#line 224 "cunmhr.f"
    if (left) {
#line 225 "cunmhr.f"
	nq = *m;
#line 226 "cunmhr.f"
	nw = *n;
#line 227 "cunmhr.f"
    } else {
#line 228 "cunmhr.f"
	nq = *n;
#line 229 "cunmhr.f"
	nw = *m;
#line 230 "cunmhr.f"
    }
#line 231 "cunmhr.f"
    if (! left && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 232 "cunmhr.f"
	*info = -1;
#line 233 "cunmhr.f"
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "C", (ftnlen)1, (ftnlen)1)) {
#line 235 "cunmhr.f"
	*info = -2;
#line 236 "cunmhr.f"
    } else if (*m < 0) {
#line 237 "cunmhr.f"
	*info = -3;
#line 238 "cunmhr.f"
    } else if (*n < 0) {
#line 239 "cunmhr.f"
	*info = -4;
#line 240 "cunmhr.f"
    } else if (*ilo < 1 || *ilo > max(1,nq)) {
#line 241 "cunmhr.f"
	*info = -5;
#line 242 "cunmhr.f"
    } else if (*ihi < min(*ilo,nq) || *ihi > nq) {
#line 243 "cunmhr.f"
	*info = -6;
#line 244 "cunmhr.f"
    } else if (*lda < max(1,nq)) {
#line 245 "cunmhr.f"
	*info = -8;
#line 246 "cunmhr.f"
    } else if (*ldc < max(1,*m)) {
#line 247 "cunmhr.f"
	*info = -11;
#line 248 "cunmhr.f"
    } else if (*lwork < max(1,nw) && ! lquery) {
#line 249 "cunmhr.f"
	*info = -13;
#line 250 "cunmhr.f"
    }

#line 252 "cunmhr.f"
    if (*info == 0) {
#line 253 "cunmhr.f"
	if (left) {
/* Writing concatenation */
#line 254 "cunmhr.f"
	    i__1[0] = 1, a__1[0] = side;
#line 254 "cunmhr.f"
	    i__1[1] = 1, a__1[1] = trans;
#line 254 "cunmhr.f"
	    s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
#line 254 "cunmhr.f"
	    nb = ilaenv_(&c__1, "CUNMQR", ch__1, &nh, n, &nh, &c_n1, (ftnlen)
		    6, (ftnlen)2);
#line 255 "cunmhr.f"
	} else {
/* Writing concatenation */
#line 256 "cunmhr.f"
	    i__1[0] = 1, a__1[0] = side;
#line 256 "cunmhr.f"
	    i__1[1] = 1, a__1[1] = trans;
#line 256 "cunmhr.f"
	    s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
#line 256 "cunmhr.f"
	    nb = ilaenv_(&c__1, "CUNMQR", ch__1, m, &nh, &nh, &c_n1, (ftnlen)
		    6, (ftnlen)2);
#line 257 "cunmhr.f"
	}
#line 258 "cunmhr.f"
	lwkopt = max(1,nw) * nb;
#line 259 "cunmhr.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 260 "cunmhr.f"
    }

#line 262 "cunmhr.f"
    if (*info != 0) {
#line 263 "cunmhr.f"
	i__2 = -(*info);
#line 263 "cunmhr.f"
	xerbla_("CUNMHR", &i__2, (ftnlen)6);
#line 264 "cunmhr.f"
	return 0;
#line 265 "cunmhr.f"
    } else if (lquery) {
#line 266 "cunmhr.f"
	return 0;
#line 267 "cunmhr.f"
    }

/*     Quick return if possible */

#line 271 "cunmhr.f"
    if (*m == 0 || *n == 0 || nh == 0) {
#line 272 "cunmhr.f"
	work[1].r = 1., work[1].i = 0.;
#line 273 "cunmhr.f"
	return 0;
#line 274 "cunmhr.f"
    }

#line 276 "cunmhr.f"
    if (left) {
#line 277 "cunmhr.f"
	mi = nh;
#line 278 "cunmhr.f"
	ni = *n;
#line 279 "cunmhr.f"
	i1 = *ilo + 1;
#line 280 "cunmhr.f"
	i2 = 1;
#line 281 "cunmhr.f"
    } else {
#line 282 "cunmhr.f"
	mi = *m;
#line 283 "cunmhr.f"
	ni = nh;
#line 284 "cunmhr.f"
	i1 = 1;
#line 285 "cunmhr.f"
	i2 = *ilo + 1;
#line 286 "cunmhr.f"
    }

#line 288 "cunmhr.f"
    cunmqr_(side, trans, &mi, &ni, &nh, &a[*ilo + 1 + *ilo * a_dim1], lda, &
	    tau[*ilo], &c__[i1 + i2 * c_dim1], ldc, &work[1], lwork, &iinfo, (
	    ftnlen)1, (ftnlen)1);

#line 291 "cunmhr.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 292 "cunmhr.f"
    return 0;

/*     End of CUNMHR */

} /* cunmhr_ */

