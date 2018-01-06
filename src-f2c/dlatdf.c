#line 1 "dlatdf.f"
/* dlatdf.f -- translated by f2c (version 20100827).
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

#line 1 "dlatdf.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b23 = 1.;
static doublereal c_b37 = -1.;

/* > \brief \b DLATDF uses the LU factorization of the n-by-n matrix computed by sgetc2 and computes a contrib
ution to the reciprocal Dif-estimate. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLATDF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlatdf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlatdf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlatdf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLATDF( IJOB, N, Z, LDZ, RHS, RDSUM, RDSCAL, IPIV, */
/*                          JPIV ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            IJOB, LDZ, N */
/*       DOUBLE PRECISION   RDSCAL, RDSUM */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ), JPIV( * ) */
/*       DOUBLE PRECISION   RHS( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLATDF uses the LU factorization of the n-by-n matrix Z computed by */
/* > DGETC2 and computes a contribution to the reciprocal Dif-estimate */
/* > by solving Z * x = b for x, and choosing the r.h.s. b such that */
/* > the norm of x is as large as possible. On entry RHS = b holds the */
/* > contribution from earlier solved sub-systems, and on return RHS = x. */
/* > */
/* > The factorization of Z returned by DGETC2 has the form Z = P*L*U*Q, */
/* > where P and Q are permutation matrices. L is lower triangular with */
/* > unit diagonal elements and U is upper triangular. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] IJOB */
/* > \verbatim */
/* >          IJOB is INTEGER */
/* >          IJOB = 2: First compute an approximative null-vector e */
/* >              of Z using DGECON, e is normalized and solve for */
/* >              Zx = +-e - f with the sign giving the greater value */
/* >              of 2-norm(x). About 5 times as expensive as Default. */
/* >          IJOB .ne. 2: Local look ahead strategy where all entries of */
/* >              the r.h.s. b is chosen as either +1 or -1 (Default). */
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
/* >          Z is DOUBLE PRECISION array, dimension (LDZ, N) */
/* >          On entry, the LU part of the factorization of the n-by-n */
/* >          matrix Z computed by DGETC2:  Z = P * L * U * Q */
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
/* >          RHS is DOUBLE PRECISION array, dimension (N) */
/* >          On entry, RHS contains contributions from other subsystems. */
/* >          On exit, RHS contains the solution of the subsystem with */
/* >          entries acoording to the value of IJOB (see above). */
/* > \endverbatim */
/* > */
/* > \param[in,out] RDSUM */
/* > \verbatim */
/* >          RDSUM is DOUBLE PRECISION */
/* >          On entry, the sum of squares of computed contributions to */
/* >          the Dif-estimate under computation by DTGSYL, where the */
/* >          scaling factor RDSCAL (see below) has been factored out. */
/* >          On exit, the corresponding sum of squares updated with the */
/* >          contributions from the current sub-system. */
/* >          If TRANS = 'T' RDSUM is not touched. */
/* >          NOTE: RDSUM only makes sense when DTGSY2 is called by STGSYL. */
/* > \endverbatim */
/* > */
/* > \param[in,out] RDSCAL */
/* > \verbatim */
/* >          RDSCAL is DOUBLE PRECISION */
/* >          On entry, scaling factor used to prevent overflow in RDSUM. */
/* >          On exit, RDSCAL is updated w.r.t. the current contributions */
/* >          in RDSUM. */
/* >          If TRANS = 'T', RDSCAL is not touched. */
/* >          NOTE: RDSCAL only makes sense when DTGSY2 is called by */
/* >                DTGSYL. */
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

/* > \ingroup doubleOTHERauxiliary */

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
/* > \verbatim */
/* > */
/* > */
/* >  [1] Bo Kagstrom and Lars Westin, */
/* >      Generalized Schur Methods with Condition Estimators for */
/* >      Solving the Generalized Sylvester Equation, IEEE Transactions */
/* >      on Automatic Control, Vol. 34, No. 7, July 1989, pp 745-751. */
/* > */
/* >  [2] Peter Poromaa, */
/* >      On Efficient and Robust Estimators for the Separation */
/* >      between two Regular Matrix Pairs with Applications in */
/* >      Condition Estimation. Report IMINF-95.05, Departement of */
/* >      Computing Science, Umea University, S-901 87 Umea, Sweden, 1995. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dlatdf_(integer *ijob, integer *n, doublereal *z__, 
	integer *ldz, doublereal *rhs, doublereal *rdsum, doublereal *rdscal, 
	integer *ipiv, integer *jpiv)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k;
    static doublereal bm, bp, xm[8], xp[8];
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer info;
    static doublereal temp, work[32];
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern doublereal dasum_(integer *, doublereal *, integer *);
    static doublereal pmone;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static doublereal sminu;
    static integer iwork[8];
    static doublereal splus;
    extern /* Subroutine */ int dgesc2_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, doublereal *), dgecon_(char *,
	     integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, ftnlen), dlassq_(integer *, 
	    doublereal *, integer *, doublereal *, doublereal *), dlaswp_(
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, integer *);


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

#line 217 "dlatdf.f"
    /* Parameter adjustments */
#line 217 "dlatdf.f"
    z_dim1 = *ldz;
#line 217 "dlatdf.f"
    z_offset = 1 + z_dim1;
#line 217 "dlatdf.f"
    z__ -= z_offset;
#line 217 "dlatdf.f"
    --rhs;
#line 217 "dlatdf.f"
    --ipiv;
#line 217 "dlatdf.f"
    --jpiv;
#line 217 "dlatdf.f"

#line 217 "dlatdf.f"
    /* Function Body */
#line 217 "dlatdf.f"
    if (*ijob != 2) {

/*        Apply permutations IPIV to RHS */

#line 221 "dlatdf.f"
	i__1 = *n - 1;
#line 221 "dlatdf.f"
	dlaswp_(&c__1, &rhs[1], ldz, &c__1, &i__1, &ipiv[1], &c__1);

/*        Solve for L-part choosing RHS either to +1 or -1. */

#line 225 "dlatdf.f"
	pmone = -1.;

#line 227 "dlatdf.f"
	i__1 = *n - 1;
#line 227 "dlatdf.f"
	for (j = 1; j <= i__1; ++j) {
#line 228 "dlatdf.f"
	    bp = rhs[j] + 1.;
#line 229 "dlatdf.f"
	    bm = rhs[j] - 1.;
#line 230 "dlatdf.f"
	    splus = 1.;

/*           Look-ahead for L-part RHS(1:N-1) = + or -1, SPLUS and */
/*           SMIN computed more efficiently than in BSOLVE [1]. */

#line 235 "dlatdf.f"
	    i__2 = *n - j;
#line 235 "dlatdf.f"
	    splus += ddot_(&i__2, &z__[j + 1 + j * z_dim1], &c__1, &z__[j + 1 
		    + j * z_dim1], &c__1);
#line 236 "dlatdf.f"
	    i__2 = *n - j;
#line 236 "dlatdf.f"
	    sminu = ddot_(&i__2, &z__[j + 1 + j * z_dim1], &c__1, &rhs[j + 1],
		     &c__1);
#line 237 "dlatdf.f"
	    splus *= rhs[j];
#line 238 "dlatdf.f"
	    if (splus > sminu) {
#line 239 "dlatdf.f"
		rhs[j] = bp;
#line 240 "dlatdf.f"
	    } else if (sminu > splus) {
#line 241 "dlatdf.f"
		rhs[j] = bm;
#line 242 "dlatdf.f"
	    } else {

/*              In this case the updating sums are equal and we can */
/*              choose RHS(J) +1 or -1. The first time this happens */
/*              we choose -1, thereafter +1. This is a simple way to */
/*              get good estimates of matrices like Byers well-known */
/*              example (see [1]). (Not done in BSOLVE.) */

#line 250 "dlatdf.f"
		rhs[j] += pmone;
#line 251 "dlatdf.f"
		pmone = 1.;
#line 252 "dlatdf.f"
	    }

/*           Compute the remaining r.h.s. */

#line 256 "dlatdf.f"
	    temp = -rhs[j];
#line 257 "dlatdf.f"
	    i__2 = *n - j;
#line 257 "dlatdf.f"
	    daxpy_(&i__2, &temp, &z__[j + 1 + j * z_dim1], &c__1, &rhs[j + 1],
		     &c__1);

#line 259 "dlatdf.f"
/* L10: */
#line 259 "dlatdf.f"
	}

/*        Solve for U-part, look-ahead for RHS(N) = +-1. This is not done */
/*        in BSOLVE and will hopefully give us a better estimate because */
/*        any ill-conditioning of the original matrix is transfered to U */
/*        and not to L. U(N, N) is an approximation to sigma_min(LU). */

#line 266 "dlatdf.f"
	i__1 = *n - 1;
#line 266 "dlatdf.f"
	dcopy_(&i__1, &rhs[1], &c__1, xp, &c__1);
#line 267 "dlatdf.f"
	xp[*n - 1] = rhs[*n] + 1.;
#line 268 "dlatdf.f"
	rhs[*n] += -1.;
#line 269 "dlatdf.f"
	splus = 0.;
#line 270 "dlatdf.f"
	sminu = 0.;
#line 271 "dlatdf.f"
	for (i__ = *n; i__ >= 1; --i__) {
#line 272 "dlatdf.f"
	    temp = 1. / z__[i__ + i__ * z_dim1];
#line 273 "dlatdf.f"
	    xp[i__ - 1] *= temp;
#line 274 "dlatdf.f"
	    rhs[i__] *= temp;
#line 275 "dlatdf.f"
	    i__1 = *n;
#line 275 "dlatdf.f"
	    for (k = i__ + 1; k <= i__1; ++k) {
#line 276 "dlatdf.f"
		xp[i__ - 1] -= xp[k - 1] * (z__[i__ + k * z_dim1] * temp);
#line 277 "dlatdf.f"
		rhs[i__] -= rhs[k] * (z__[i__ + k * z_dim1] * temp);
#line 278 "dlatdf.f"
/* L20: */
#line 278 "dlatdf.f"
	    }
#line 279 "dlatdf.f"
	    splus += (d__1 = xp[i__ - 1], abs(d__1));
#line 280 "dlatdf.f"
	    sminu += (d__1 = rhs[i__], abs(d__1));
#line 281 "dlatdf.f"
/* L30: */
#line 281 "dlatdf.f"
	}
#line 282 "dlatdf.f"
	if (splus > sminu) {
#line 282 "dlatdf.f"
	    dcopy_(n, xp, &c__1, &rhs[1], &c__1);
#line 282 "dlatdf.f"
	}

/*        Apply the permutations JPIV to the computed solution (RHS) */

#line 287 "dlatdf.f"
	i__1 = *n - 1;
#line 287 "dlatdf.f"
	dlaswp_(&c__1, &rhs[1], ldz, &c__1, &i__1, &jpiv[1], &c_n1);

/*        Compute the sum of squares */

#line 291 "dlatdf.f"
	dlassq_(n, &rhs[1], &c__1, rdscal, rdsum);

#line 293 "dlatdf.f"
    } else {

/*        IJOB = 2, Compute approximate nullvector XM of Z */

#line 297 "dlatdf.f"
	dgecon_("I", n, &z__[z_offset], ldz, &c_b23, &temp, work, iwork, &
		info, (ftnlen)1);
#line 298 "dlatdf.f"
	dcopy_(n, &work[*n], &c__1, xm, &c__1);

/*        Compute RHS */

#line 302 "dlatdf.f"
	i__1 = *n - 1;
#line 302 "dlatdf.f"
	dlaswp_(&c__1, xm, ldz, &c__1, &i__1, &ipiv[1], &c_n1);
#line 303 "dlatdf.f"
	temp = 1. / sqrt(ddot_(n, xm, &c__1, xm, &c__1));
#line 304 "dlatdf.f"
	dscal_(n, &temp, xm, &c__1);
#line 305 "dlatdf.f"
	dcopy_(n, xm, &c__1, xp, &c__1);
#line 306 "dlatdf.f"
	daxpy_(n, &c_b23, &rhs[1], &c__1, xp, &c__1);
#line 307 "dlatdf.f"
	daxpy_(n, &c_b37, xm, &c__1, &rhs[1], &c__1);
#line 308 "dlatdf.f"
	dgesc2_(n, &z__[z_offset], ldz, &rhs[1], &ipiv[1], &jpiv[1], &temp);
#line 309 "dlatdf.f"
	dgesc2_(n, &z__[z_offset], ldz, xp, &ipiv[1], &jpiv[1], &temp);
#line 310 "dlatdf.f"
	if (dasum_(n, xp, &c__1) > dasum_(n, &rhs[1], &c__1)) {
#line 310 "dlatdf.f"
	    dcopy_(n, xp, &c__1, &rhs[1], &c__1);
#line 310 "dlatdf.f"
	}

/*        Compute the sum of squares */

#line 315 "dlatdf.f"
	dlassq_(n, &rhs[1], &c__1, rdscal, rdsum);

#line 317 "dlatdf.f"
    }

#line 319 "dlatdf.f"
    return 0;

/*     End of DLATDF */

} /* dlatdf_ */

