#line 1 "dlasq1.f"
/* dlasq1.f -- translated by f2c (version 20100827).
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

#line 1 "dlasq1.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;
static integer c__0 = 0;

/* > \brief \b DLASQ1 computes the singular values of a real square bidiagonal matrix. Used by sbdsqr. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLASQ1 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasq1.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasq1.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasq1.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLASQ1( N, D, E, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   D( * ), E( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLASQ1 computes the singular values of a real N-by-N bidiagonal */
/* > matrix with diagonal D and off-diagonal E. The singular values */
/* > are computed to high relative accuracy, in the absence of */
/* > denormalization, underflow and overflow. The algorithm was first */
/* > presented in */
/* > */
/* > "Accurate singular values and differential qd algorithms" by K. V. */
/* > Fernando and B. N. Parlett, Numer. Math., Vol-67, No. 2, pp. 191-230, */
/* > 1994, */
/* > */
/* > and the present implementation is described in "An implementation of */
/* > the dqds Algorithm (Positive Case)", LAPACK Working Note. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >        The number of rows and columns in the matrix. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >        On entry, D contains the diagonal elements of the */
/* >        bidiagonal matrix whose SVD is desired. On normal exit, */
/* >        D contains the singular values in decreasing order. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array, dimension (N) */
/* >        On entry, elements E(1:N-1) contain the off-diagonal elements */
/* >        of the bidiagonal matrix whose SVD is desired. */
/* >        On exit, E is overwritten. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (4*N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >        = 0: successful exit */
/* >        < 0: if INFO = -i, the i-th argument had an illegal value */
/* >        > 0: the algorithm failed */
/* >             = 1, a split was marked by a positive value in E */
/* >             = 2, current block of Z not diagonalized after 100*N */
/* >                  iterations (in inner while loop)  On exit D and E */
/* >                  represent a matrix with the same singular values */
/* >                  which the calling subroutine could use to finish the */
/* >                  computation, or even feed back into DLASQ1 */
/* >             = 3, termination criterion of outer while loop not met */
/* >                  (program created more than N unreduced blocks) */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup auxOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int dlasq1_(integer *n, doublereal *d__, doublereal *e, 
	doublereal *work, integer *info)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal eps;
    extern /* Subroutine */ int dlas2_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *);
    static doublereal scale;
    static integer iinfo;
    static doublereal sigmn;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal sigmx;
    extern /* Subroutine */ int dlasq2_(integer *, doublereal *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), dlasrt_(
	    char *, integer *, doublereal *, integer *, ftnlen);


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
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 145 "dlasq1.f"
    /* Parameter adjustments */
#line 145 "dlasq1.f"
    --work;
#line 145 "dlasq1.f"
    --e;
#line 145 "dlasq1.f"
    --d__;
#line 145 "dlasq1.f"

#line 145 "dlasq1.f"
    /* Function Body */
#line 145 "dlasq1.f"
    *info = 0;
#line 146 "dlasq1.f"
    if (*n < 0) {
#line 147 "dlasq1.f"
	*info = -2;
#line 148 "dlasq1.f"
	i__1 = -(*info);
#line 148 "dlasq1.f"
	xerbla_("DLASQ1", &i__1, (ftnlen)6);
#line 149 "dlasq1.f"
	return 0;
#line 150 "dlasq1.f"
    } else if (*n == 0) {
#line 151 "dlasq1.f"
	return 0;
#line 152 "dlasq1.f"
    } else if (*n == 1) {
#line 153 "dlasq1.f"
	d__[1] = abs(d__[1]);
#line 154 "dlasq1.f"
	return 0;
#line 155 "dlasq1.f"
    } else if (*n == 2) {
#line 156 "dlasq1.f"
	dlas2_(&d__[1], &e[1], &d__[2], &sigmn, &sigmx);
#line 157 "dlasq1.f"
	d__[1] = sigmx;
#line 158 "dlasq1.f"
	d__[2] = sigmn;
#line 159 "dlasq1.f"
	return 0;
#line 160 "dlasq1.f"
    }

/*     Estimate the largest singular value. */

#line 164 "dlasq1.f"
    sigmx = 0.;
#line 165 "dlasq1.f"
    i__1 = *n - 1;
#line 165 "dlasq1.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 166 "dlasq1.f"
	d__[i__] = (d__1 = d__[i__], abs(d__1));
/* Computing MAX */
#line 167 "dlasq1.f"
	d__2 = sigmx, d__3 = (d__1 = e[i__], abs(d__1));
#line 167 "dlasq1.f"
	sigmx = max(d__2,d__3);
#line 168 "dlasq1.f"
/* L10: */
#line 168 "dlasq1.f"
    }
#line 169 "dlasq1.f"
    d__[*n] = (d__1 = d__[*n], abs(d__1));

/*     Early return if SIGMX is zero (matrix is already diagonal). */

#line 173 "dlasq1.f"
    if (sigmx == 0.) {
#line 174 "dlasq1.f"
	dlasrt_("D", n, &d__[1], &iinfo, (ftnlen)1);
#line 175 "dlasq1.f"
	return 0;
#line 176 "dlasq1.f"
    }

#line 178 "dlasq1.f"
    i__1 = *n;
#line 178 "dlasq1.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 179 "dlasq1.f"
	d__1 = sigmx, d__2 = d__[i__];
#line 179 "dlasq1.f"
	sigmx = max(d__1,d__2);
#line 180 "dlasq1.f"
/* L20: */
#line 180 "dlasq1.f"
    }

/*     Copy D and E into WORK (in the Z format) and scale (squaring the */
/*     input data makes scaling by a power of the radix pointless). */

#line 185 "dlasq1.f"
    eps = dlamch_("Precision", (ftnlen)9);
#line 186 "dlasq1.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 187 "dlasq1.f"
    scale = sqrt(eps / safmin);
#line 188 "dlasq1.f"
    dcopy_(n, &d__[1], &c__1, &work[1], &c__2);
#line 189 "dlasq1.f"
    i__1 = *n - 1;
#line 189 "dlasq1.f"
    dcopy_(&i__1, &e[1], &c__1, &work[2], &c__2);
#line 190 "dlasq1.f"
    i__1 = (*n << 1) - 1;
#line 190 "dlasq1.f"
    i__2 = (*n << 1) - 1;
#line 190 "dlasq1.f"
    dlascl_("G", &c__0, &c__0, &sigmx, &scale, &i__1, &c__1, &work[1], &i__2, 
	    &iinfo, (ftnlen)1);

/*     Compute the q's and e's. */

#line 195 "dlasq1.f"
    i__1 = (*n << 1) - 1;
#line 195 "dlasq1.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
#line 196 "dlasq1.f"
	d__1 = work[i__];
#line 196 "dlasq1.f"
	work[i__] = d__1 * d__1;
#line 197 "dlasq1.f"
/* L30: */
#line 197 "dlasq1.f"
    }
#line 198 "dlasq1.f"
    work[*n * 2] = 0.;

#line 200 "dlasq1.f"
    dlasq2_(n, &work[1], info);

#line 202 "dlasq1.f"
    if (*info == 0) {
#line 203 "dlasq1.f"
	i__1 = *n;
#line 203 "dlasq1.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 204 "dlasq1.f"
	    d__[i__] = sqrt(work[i__]);
#line 205 "dlasq1.f"
/* L40: */
#line 205 "dlasq1.f"
	}
#line 206 "dlasq1.f"
	dlascl_("G", &c__0, &c__0, &scale, &sigmx, n, &c__1, &d__[1], n, &
		iinfo, (ftnlen)1);
#line 207 "dlasq1.f"
    } else if (*info == 2) {

/*     Maximum number of iterations exceeded.  Move data from WORK */
/*     into D and E so the calling subroutine can try to finish */

#line 212 "dlasq1.f"
	i__1 = *n;
#line 212 "dlasq1.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 213 "dlasq1.f"
	    d__[i__] = sqrt(work[(i__ << 1) - 1]);
#line 214 "dlasq1.f"
	    e[i__] = sqrt(work[i__ * 2]);
#line 215 "dlasq1.f"
	}
#line 216 "dlasq1.f"
	dlascl_("G", &c__0, &c__0, &scale, &sigmx, n, &c__1, &d__[1], n, &
		iinfo, (ftnlen)1);
#line 217 "dlasq1.f"
	dlascl_("G", &c__0, &c__0, &scale, &sigmx, n, &c__1, &e[1], n, &iinfo,
		 (ftnlen)1);
#line 218 "dlasq1.f"
    }

#line 220 "dlasq1.f"
    return 0;

/*     End of DLASQ1 */

} /* dlasq1_ */

