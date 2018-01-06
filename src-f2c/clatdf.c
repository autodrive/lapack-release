#line 1 "clatdf.f"
/* clatdf.f -- translated by f2c (version 20100827).
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

#line 1 "clatdf.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b24 = 1.;

/* > \brief \b CLATDF uses the LU factorization of the n-by-n matrix computed by sgetc2 and computes a contrib
ution to the reciprocal Dif-estimate. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLATDF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clatdf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clatdf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clatdf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLATDF( IJOB, N, Z, LDZ, RHS, RDSUM, RDSCAL, IPIV, */
/*                          JPIV ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            IJOB, LDZ, N */
/*       REAL               RDSCAL, RDSUM */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ), JPIV( * ) */
/*       COMPLEX            RHS( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLATDF computes the contribution to the reciprocal Dif-estimate */
/* > by solving for x in Z * x = b, where b is chosen such that the norm */
/* > of x is as large as possible. It is assumed that LU decomposition */
/* > of Z has been computed by CGETC2. On entry RHS = f holds the */
/* > contribution from earlier solved sub-systems, and on return RHS = x. */
/* > */
/* > The factorization of Z returned by CGETC2 has the form */
/* > Z = P * L * U * Q, where P and Q are permutation matrices. L is lower */
/* > triangular with unit diagonal elements and U is upper triangular. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] IJOB */
/* > \verbatim */
/* >          IJOB is INTEGER */
/* >          IJOB = 2: First compute an approximative null-vector e */
/* >              of Z using CGECON, e is normalized and solve for */
/* >              Zx = +-e - f with the sign giving the greater value of */
/* >              2-norm(x).  About 5 times as expensive as Default. */
/* >          IJOB .ne. 2: Local look ahead strategy where */
/* >              all entries of the r.h.s. b is chosen as either +1 or */
/* >              -1.  Default. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix Z. */
/* > \endverbatim */
/* > */
/* > \param[in] Z */
/* > \verbatim */
/* >          Z is COMPLEX array, dimension (LDZ, N) */
/* >          On entry, the LU part of the factorization of the n-by-n */
/* >          matrix Z computed by CGETC2:  Z = P * L * U * Q */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* >          LDZ is INTEGER */
/* >          The leading dimension of the array Z.  LDA >= max(1, N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] RHS */
/* > \verbatim */
/* >          RHS is COMPLEX array, dimension (N). */
/* >          On entry, RHS contains contributions from other subsystems. */
/* >          On exit, RHS contains the solution of the subsystem with */
/* >          entries according to the value of IJOB (see above). */
/* > \endverbatim */
/* > */
/* > \param[in,out] RDSUM */
/* > \verbatim */
/* >          RDSUM is REAL */
/* >          On entry, the sum of squares of computed contributions to */
/* >          the Dif-estimate under computation by CTGSYL, where the */
/* >          scaling factor RDSCAL (see below) has been factored out. */
/* >          On exit, the corresponding sum of squares updated with the */
/* >          contributions from the current sub-system. */
/* >          If TRANS = 'T' RDSUM is not touched. */
/* >          NOTE: RDSUM only makes sense when CTGSY2 is called by CTGSYL. */
/* > \endverbatim */
/* > */
/* > \param[in,out] RDSCAL */
/* > \verbatim */
/* >          RDSCAL is REAL */
/* >          On entry, scaling factor used to prevent overflow in RDSUM. */
/* >          On exit, RDSCAL is updated w.r.t. the current contributions */
/* >          in RDSUM. */
/* >          If TRANS = 'T', RDSCAL is not touched. */
/* >          NOTE: RDSCAL only makes sense when CTGSY2 is called by */
/* >          CTGSYL. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N). */
/* >          The pivot indices; for 1 <= i <= N, row i of the */
/* >          matrix has been interchanged with row IPIV(i). */
/* > \endverbatim */
/* > */
/* > \param[in] JPIV */
/* > \verbatim */
/* >          JPIV is INTEGER array, dimension (N). */
/* >          The pivot indices; for 1 <= j <= N, column j of the */
/* >          matrix has been interchanged with column JPIV(j). */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2016 */

/* > \ingroup complexOTHERauxiliary */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* >  This routine is a further developed implementation of algorithm */
/* >  BSOLVE in [1] using complete pivoting in the LU factorization. */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Bo Kagstrom and Peter Poromaa, Department of Computing Science, */
/* >     Umea University, S-901 87 Umea, Sweden. */

/* > \par References: */
/*  ================ */
/* > */
/* >   [1]   Bo Kagstrom and Lars Westin, */
/* >         Generalized Schur Methods with Condition Estimators for */
/* >         Solving the Generalized Sylvester Equation, IEEE Transactions */
/* >         on Automatic Control, Vol. 34, No. 7, July 1989, pp 745-751. */
/* > */
/* >   [2]   Peter Poromaa, */
/* >         On Efficient and Robust Estimators for the Separation */
/* >         between two Regular Matrix Pairs with Applications in */
/* >         Condition Estimation. Report UMINF-95.05, Department of */
/* >         Computing Science, Umea University, S-901 87 Umea, Sweden, */
/* >         1995. */

/*  ===================================================================== */
/* Subroutine */ int clatdf_(integer *ijob, integer *n, doublecomplex *z__, 
	integer *ldz, doublecomplex *rhs, doublereal *rdsum, doublereal *
	rdscal, integer *ipiv, integer *jpiv)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);
    double z_abs(doublecomplex *);
    void z_sqrt(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, k;
    static doublecomplex bm, bp, xm[2], xp[2];
    static integer info;
    static doublecomplex temp, work[8];
    extern /* Subroutine */ int cscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *);
    static doublereal scale;
    extern /* Double Complex */ VOID cdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int ccopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static doublecomplex pmone;
    extern /* Subroutine */ int caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static doublereal rtemp, sminu, rwork[2], splus;
    extern /* Subroutine */ int cgesc2_(integer *, doublecomplex *, integer *,
	     doublecomplex *, integer *, integer *, doublereal *), cgecon_(
	    char *, integer *, doublecomplex *, integer *, doublereal *, 
	    doublereal *, doublecomplex *, doublereal *, integer *, ftnlen), 
	    classq_(integer *, doublecomplex *, integer *, doublereal *, 
	    doublereal *), claswp_(integer *, doublecomplex *, integer *, 
	    integer *, integer *, integer *, integer *);
    extern doublereal scasum_(integer *, doublecomplex *, integer *);


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 219 "clatdf.f"
    /* Parameter adjustments */
#line 219 "clatdf.f"
    z_dim1 = *ldz;
#line 219 "clatdf.f"
    z_offset = 1 + z_dim1;
#line 219 "clatdf.f"
    z__ -= z_offset;
#line 219 "clatdf.f"
    --rhs;
#line 219 "clatdf.f"
    --ipiv;
#line 219 "clatdf.f"
    --jpiv;
#line 219 "clatdf.f"

#line 219 "clatdf.f"
    /* Function Body */
#line 219 "clatdf.f"
    if (*ijob != 2) {

/*        Apply permutations IPIV to RHS */

#line 223 "clatdf.f"
	i__1 = *n - 1;
#line 223 "clatdf.f"
	claswp_(&c__1, &rhs[1], ldz, &c__1, &i__1, &ipiv[1], &c__1);

/*        Solve for L-part choosing RHS either to +1 or -1. */

#line 227 "clatdf.f"
	z__1.r = -1., z__1.i = -0.;
#line 227 "clatdf.f"
	pmone.r = z__1.r, pmone.i = z__1.i;
#line 228 "clatdf.f"
	i__1 = *n - 1;
#line 228 "clatdf.f"
	for (j = 1; j <= i__1; ++j) {
#line 229 "clatdf.f"
	    i__2 = j;
#line 229 "clatdf.f"
	    z__1.r = rhs[i__2].r + 1., z__1.i = rhs[i__2].i + 0.;
#line 229 "clatdf.f"
	    bp.r = z__1.r, bp.i = z__1.i;
#line 230 "clatdf.f"
	    i__2 = j;
#line 230 "clatdf.f"
	    z__1.r = rhs[i__2].r - 1., z__1.i = rhs[i__2].i - 0.;
#line 230 "clatdf.f"
	    bm.r = z__1.r, bm.i = z__1.i;
#line 231 "clatdf.f"
	    splus = 1.;

/*           Lockahead for L- part RHS(1:N-1) = +-1 */
/*           SPLUS and SMIN computed more efficiently than in BSOLVE[1]. */

#line 236 "clatdf.f"
	    i__2 = *n - j;
#line 236 "clatdf.f"
	    cdotc_(&z__1, &i__2, &z__[j + 1 + j * z_dim1], &c__1, &z__[j + 1 
		    + j * z_dim1], &c__1);
#line 236 "clatdf.f"
	    splus += z__1.r;
#line 238 "clatdf.f"
	    i__2 = *n - j;
#line 238 "clatdf.f"
	    cdotc_(&z__1, &i__2, &z__[j + 1 + j * z_dim1], &c__1, &rhs[j + 1],
		     &c__1);
#line 238 "clatdf.f"
	    sminu = z__1.r;
#line 239 "clatdf.f"
	    i__2 = j;
#line 239 "clatdf.f"
	    splus *= rhs[i__2].r;
#line 240 "clatdf.f"
	    if (splus > sminu) {
#line 241 "clatdf.f"
		i__2 = j;
#line 241 "clatdf.f"
		rhs[i__2].r = bp.r, rhs[i__2].i = bp.i;
#line 242 "clatdf.f"
	    } else if (sminu > splus) {
#line 243 "clatdf.f"
		i__2 = j;
#line 243 "clatdf.f"
		rhs[i__2].r = bm.r, rhs[i__2].i = bm.i;
#line 244 "clatdf.f"
	    } else {

/*              In this case the updating sums are equal and we can */
/*              choose RHS(J) +1 or -1. The first time this happens we */
/*              choose -1, thereafter +1. This is a simple way to get */
/*              good estimates of matrices like Byers well-known example */
/*              (see [1]). (Not done in BSOLVE.) */

#line 252 "clatdf.f"
		i__2 = j;
#line 252 "clatdf.f"
		i__3 = j;
#line 252 "clatdf.f"
		z__1.r = rhs[i__3].r + pmone.r, z__1.i = rhs[i__3].i + 
			pmone.i;
#line 252 "clatdf.f"
		rhs[i__2].r = z__1.r, rhs[i__2].i = z__1.i;
#line 253 "clatdf.f"
		pmone.r = 1., pmone.i = 0.;
#line 254 "clatdf.f"
	    }

/*           Compute the remaining r.h.s. */

#line 258 "clatdf.f"
	    i__2 = j;
#line 258 "clatdf.f"
	    z__1.r = -rhs[i__2].r, z__1.i = -rhs[i__2].i;
#line 258 "clatdf.f"
	    temp.r = z__1.r, temp.i = z__1.i;
#line 259 "clatdf.f"
	    i__2 = *n - j;
#line 259 "clatdf.f"
	    caxpy_(&i__2, &temp, &z__[j + 1 + j * z_dim1], &c__1, &rhs[j + 1],
		     &c__1);
#line 260 "clatdf.f"
/* L10: */
#line 260 "clatdf.f"
	}

/*        Solve for U- part, lockahead for RHS(N) = +-1. This is not done */
/*        In BSOLVE and will hopefully give us a better estimate because */
/*        any ill-conditioning of the original matrix is transfered to U */
/*        and not to L. U(N, N) is an approximation to sigma_min(LU). */

#line 267 "clatdf.f"
	i__1 = *n - 1;
#line 267 "clatdf.f"
	ccopy_(&i__1, &rhs[1], &c__1, work, &c__1);
#line 268 "clatdf.f"
	i__1 = *n - 1;
#line 268 "clatdf.f"
	i__2 = *n;
#line 268 "clatdf.f"
	z__1.r = rhs[i__2].r + 1., z__1.i = rhs[i__2].i + 0.;
#line 268 "clatdf.f"
	work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 269 "clatdf.f"
	i__1 = *n;
#line 269 "clatdf.f"
	i__2 = *n;
#line 269 "clatdf.f"
	z__1.r = rhs[i__2].r - 1., z__1.i = rhs[i__2].i - 0.;
#line 269 "clatdf.f"
	rhs[i__1].r = z__1.r, rhs[i__1].i = z__1.i;
#line 270 "clatdf.f"
	splus = 0.;
#line 271 "clatdf.f"
	sminu = 0.;
#line 272 "clatdf.f"
	for (i__ = *n; i__ >= 1; --i__) {
#line 273 "clatdf.f"
	    z_div(&z__1, &c_b1, &z__[i__ + i__ * z_dim1]);
#line 273 "clatdf.f"
	    temp.r = z__1.r, temp.i = z__1.i;
#line 274 "clatdf.f"
	    i__1 = i__ - 1;
#line 274 "clatdf.f"
	    i__2 = i__ - 1;
#line 274 "clatdf.f"
	    z__1.r = work[i__2].r * temp.r - work[i__2].i * temp.i, z__1.i = 
		    work[i__2].r * temp.i + work[i__2].i * temp.r;
#line 274 "clatdf.f"
	    work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 275 "clatdf.f"
	    i__1 = i__;
#line 275 "clatdf.f"
	    i__2 = i__;
#line 275 "clatdf.f"
	    z__1.r = rhs[i__2].r * temp.r - rhs[i__2].i * temp.i, z__1.i = 
		    rhs[i__2].r * temp.i + rhs[i__2].i * temp.r;
#line 275 "clatdf.f"
	    rhs[i__1].r = z__1.r, rhs[i__1].i = z__1.i;
#line 276 "clatdf.f"
	    i__1 = *n;
#line 276 "clatdf.f"
	    for (k = i__ + 1; k <= i__1; ++k) {
#line 277 "clatdf.f"
		i__2 = i__ - 1;
#line 277 "clatdf.f"
		i__3 = i__ - 1;
#line 277 "clatdf.f"
		i__4 = k - 1;
#line 277 "clatdf.f"
		i__5 = i__ + k * z_dim1;
#line 277 "clatdf.f"
		z__3.r = z__[i__5].r * temp.r - z__[i__5].i * temp.i, z__3.i =
			 z__[i__5].r * temp.i + z__[i__5].i * temp.r;
#line 277 "clatdf.f"
		z__2.r = work[i__4].r * z__3.r - work[i__4].i * z__3.i, 
			z__2.i = work[i__4].r * z__3.i + work[i__4].i * 
			z__3.r;
#line 277 "clatdf.f"
		z__1.r = work[i__3].r - z__2.r, z__1.i = work[i__3].i - 
			z__2.i;
#line 277 "clatdf.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 278 "clatdf.f"
		i__2 = i__;
#line 278 "clatdf.f"
		i__3 = i__;
#line 278 "clatdf.f"
		i__4 = k;
#line 278 "clatdf.f"
		i__5 = i__ + k * z_dim1;
#line 278 "clatdf.f"
		z__3.r = z__[i__5].r * temp.r - z__[i__5].i * temp.i, z__3.i =
			 z__[i__5].r * temp.i + z__[i__5].i * temp.r;
#line 278 "clatdf.f"
		z__2.r = rhs[i__4].r * z__3.r - rhs[i__4].i * z__3.i, z__2.i =
			 rhs[i__4].r * z__3.i + rhs[i__4].i * z__3.r;
#line 278 "clatdf.f"
		z__1.r = rhs[i__3].r - z__2.r, z__1.i = rhs[i__3].i - z__2.i;
#line 278 "clatdf.f"
		rhs[i__2].r = z__1.r, rhs[i__2].i = z__1.i;
#line 279 "clatdf.f"
/* L20: */
#line 279 "clatdf.f"
	    }
#line 280 "clatdf.f"
	    splus += z_abs(&work[i__ - 1]);
#line 281 "clatdf.f"
	    sminu += z_abs(&rhs[i__]);
#line 282 "clatdf.f"
/* L30: */
#line 282 "clatdf.f"
	}
#line 283 "clatdf.f"
	if (splus > sminu) {
#line 283 "clatdf.f"
	    ccopy_(n, work, &c__1, &rhs[1], &c__1);
#line 283 "clatdf.f"
	}

/*        Apply the permutations JPIV to the computed solution (RHS) */

#line 288 "clatdf.f"
	i__1 = *n - 1;
#line 288 "clatdf.f"
	claswp_(&c__1, &rhs[1], ldz, &c__1, &i__1, &jpiv[1], &c_n1);

/*        Compute the sum of squares */

#line 292 "clatdf.f"
	classq_(n, &rhs[1], &c__1, rdscal, rdsum);
#line 293 "clatdf.f"
	return 0;
#line 294 "clatdf.f"
    }

/*     ENTRY IJOB = 2 */

/*     Compute approximate nullvector XM of Z */

#line 300 "clatdf.f"
    cgecon_("I", n, &z__[z_offset], ldz, &c_b24, &rtemp, work, rwork, &info, (
	    ftnlen)1);
#line 301 "clatdf.f"
    ccopy_(n, &work[*n], &c__1, xm, &c__1);

/*     Compute RHS */

#line 305 "clatdf.f"
    i__1 = *n - 1;
#line 305 "clatdf.f"
    claswp_(&c__1, xm, ldz, &c__1, &i__1, &ipiv[1], &c_n1);
#line 306 "clatdf.f"
    cdotc_(&z__3, n, xm, &c__1, xm, &c__1);
#line 306 "clatdf.f"
    z_sqrt(&z__2, &z__3);
#line 306 "clatdf.f"
    z_div(&z__1, &c_b1, &z__2);
#line 306 "clatdf.f"
    temp.r = z__1.r, temp.i = z__1.i;
#line 307 "clatdf.f"
    cscal_(n, &temp, xm, &c__1);
#line 308 "clatdf.f"
    ccopy_(n, xm, &c__1, xp, &c__1);
#line 309 "clatdf.f"
    caxpy_(n, &c_b1, &rhs[1], &c__1, xp, &c__1);
#line 310 "clatdf.f"
    z__1.r = -1., z__1.i = -0.;
#line 310 "clatdf.f"
    caxpy_(n, &z__1, xm, &c__1, &rhs[1], &c__1);
#line 311 "clatdf.f"
    cgesc2_(n, &z__[z_offset], ldz, &rhs[1], &ipiv[1], &jpiv[1], &scale);
#line 312 "clatdf.f"
    cgesc2_(n, &z__[z_offset], ldz, xp, &ipiv[1], &jpiv[1], &scale);
#line 313 "clatdf.f"
    if (scasum_(n, xp, &c__1) > scasum_(n, &rhs[1], &c__1)) {
#line 313 "clatdf.f"
	ccopy_(n, xp, &c__1, &rhs[1], &c__1);
#line 313 "clatdf.f"
    }

/*     Compute the sum of squares */

#line 318 "clatdf.f"
    classq_(n, &rhs[1], &c__1, rdscal, rdsum);
#line 319 "clatdf.f"
    return 0;

/*     End of CLATDF */

} /* clatdf_ */

