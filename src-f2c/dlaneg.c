#line 1 "dlaneg.f"
/* dlaneg.f -- translated by f2c (version 20100827).
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

#line 1 "dlaneg.f"
/* > \brief \b DLANEG computes the Sturm count. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLANEG + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaneg.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaneg.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaneg.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       INTEGER FUNCTION DLANEG( N, D, LLD, SIGMA, PIVMIN, R ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            N, R */
/*       DOUBLE PRECISION   PIVMIN, SIGMA */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   D( * ), LLD( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLANEG computes the Sturm count, the number of negative pivots */
/* > encountered while factoring tridiagonal T - sigma I = L D L^T. */
/* > This implementation works directly on the factors without forming */
/* > the tridiagonal matrix T.  The Sturm count is also the number of */
/* > eigenvalues of T less than sigma. */
/* > */
/* > This routine is called from DLARRB. */
/* > */
/* > The current routine does not use the PIVMIN parameter but rather */
/* > requires IEEE-754 propagation of Infinities and NaNs.  This */
/* > routine also has no input range restrictions but does require */
/* > default exception handling such that x/0 produces Inf when x is */
/* > non-zero, and Inf/Inf produces NaN.  For more information, see: */
/* > */
/* >   Marques, Riedy, and Voemel, "Benefits of IEEE-754 Features in */
/* >   Modern Symmetric Tridiagonal Eigensolvers," SIAM Journal on */
/* >   Scientific Computing, v28, n5, 2006.  DOI 10.1137/050641624 */
/* >   (Tech report version in LAWN 172 with the same title.) */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          The N diagonal elements of the diagonal matrix D. */
/* > \endverbatim */
/* > */
/* > \param[in] LLD */
/* > \verbatim */
/* >          LLD is DOUBLE PRECISION array, dimension (N-1) */
/* >          The (N-1) elements L(i)*L(i)*D(i). */
/* > \endverbatim */
/* > */
/* > \param[in] SIGMA */
/* > \verbatim */
/* >          SIGMA is DOUBLE PRECISION */
/* >          Shift amount in T - sigma I = L D L^T. */
/* > \endverbatim */
/* > */
/* > \param[in] PIVMIN */
/* > \verbatim */
/* >          PIVMIN is DOUBLE PRECISION */
/* >          The minimum pivot in the Sturm sequence.  May be used */
/* >          when zero pivots are encountered on non-IEEE-754 */
/* >          architectures. */
/* > \endverbatim */
/* > */
/* > \param[in] R */
/* > \verbatim */
/* >          R is INTEGER */
/* >          The twist index for the twisted factorization that is used */
/* >          for the negcount. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup OTHERauxiliary */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Osni Marques, LBNL/NERSC, USA \n */
/* >     Christof Voemel, University of California, Berkeley, USA \n */
/* >     Jason Riedy, University of California, Berkeley, USA \n */
/* > */
/*  ===================================================================== */
integer dlaneg_(integer *n, doublereal *d__, doublereal *lld, doublereal *
	sigma, doublereal *pivmin, integer *r__)
{
    /* System generated locals */
    integer ret_val, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer j;
    static doublereal p, t;
    static integer bj;
    static doublereal tmp;
    static integer neg1, neg2;
    static doublereal bsav, gamma, dplus;
    extern logical disnan_(doublereal *);
    static integer negcnt;
    static logical sawnan;
    static doublereal dminus;


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     Some architectures propagate Infinities and NaNs very slowly, so */
/*     the code computes counts in BLKLEN chunks.  Then a NaN can */
/*     propagate at most BLKLEN columns before being detected.  This is */
/*     not a general tuning parameter; it needs only to be just large */
/*     enough that the overhead is tiny in common cases. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */
#line 161 "dlaneg.f"
    /* Parameter adjustments */
#line 161 "dlaneg.f"
    --lld;
#line 161 "dlaneg.f"
    --d__;
#line 161 "dlaneg.f"

#line 161 "dlaneg.f"
    /* Function Body */
#line 161 "dlaneg.f"
    negcnt = 0;
/*     I) upper part: L D L^T - SIGMA I = L+ D+ L+^T */
#line 164 "dlaneg.f"
    t = -(*sigma);
#line 165 "dlaneg.f"
    i__1 = *r__ - 1;
#line 165 "dlaneg.f"
    for (bj = 1; bj <= i__1; bj += 128) {
#line 166 "dlaneg.f"
	neg1 = 0;
#line 167 "dlaneg.f"
	bsav = t;
/* Computing MIN */
#line 168 "dlaneg.f"
	i__3 = bj + 127, i__4 = *r__ - 1;
#line 168 "dlaneg.f"
	i__2 = min(i__3,i__4);
#line 168 "dlaneg.f"
	for (j = bj; j <= i__2; ++j) {
#line 169 "dlaneg.f"
	    dplus = d__[j] + t;
#line 170 "dlaneg.f"
	    if (dplus < 0.) {
#line 170 "dlaneg.f"
		++neg1;
#line 170 "dlaneg.f"
	    }
#line 171 "dlaneg.f"
	    tmp = t / dplus;
#line 172 "dlaneg.f"
	    t = tmp * lld[j] - *sigma;
#line 173 "dlaneg.f"
/* L21: */
#line 173 "dlaneg.f"
	}
#line 174 "dlaneg.f"
	sawnan = disnan_(&t);
/*     Run a slower version of the above loop if a NaN is detected. */
/*     A NaN should occur only with a zero pivot after an infinite */
/*     pivot.  In that case, substituting 1 for T/DPLUS is the */
/*     correct limit. */
#line 179 "dlaneg.f"
	if (sawnan) {
#line 180 "dlaneg.f"
	    neg1 = 0;
#line 181 "dlaneg.f"
	    t = bsav;
/* Computing MIN */
#line 182 "dlaneg.f"
	    i__3 = bj + 127, i__4 = *r__ - 1;
#line 182 "dlaneg.f"
	    i__2 = min(i__3,i__4);
#line 182 "dlaneg.f"
	    for (j = bj; j <= i__2; ++j) {
#line 183 "dlaneg.f"
		dplus = d__[j] + t;
#line 184 "dlaneg.f"
		if (dplus < 0.) {
#line 184 "dlaneg.f"
		    ++neg1;
#line 184 "dlaneg.f"
		}
#line 185 "dlaneg.f"
		tmp = t / dplus;
#line 186 "dlaneg.f"
		if (disnan_(&tmp)) {
#line 186 "dlaneg.f"
		    tmp = 1.;
#line 186 "dlaneg.f"
		}
#line 187 "dlaneg.f"
		t = tmp * lld[j] - *sigma;
#line 188 "dlaneg.f"
/* L22: */
#line 188 "dlaneg.f"
	    }
#line 189 "dlaneg.f"
	}
#line 190 "dlaneg.f"
	negcnt += neg1;
#line 191 "dlaneg.f"
/* L210: */
#line 191 "dlaneg.f"
    }

/*     II) lower part: L D L^T - SIGMA I = U- D- U-^T */
#line 194 "dlaneg.f"
    p = d__[*n] - *sigma;
#line 195 "dlaneg.f"
    i__1 = *r__;
#line 195 "dlaneg.f"
    for (bj = *n - 1; bj >= i__1; bj += -128) {
#line 196 "dlaneg.f"
	neg2 = 0;
#line 197 "dlaneg.f"
	bsav = p;
/* Computing MAX */
#line 198 "dlaneg.f"
	i__3 = bj - 127;
#line 198 "dlaneg.f"
	i__2 = max(i__3,*r__);
#line 198 "dlaneg.f"
	for (j = bj; j >= i__2; --j) {
#line 199 "dlaneg.f"
	    dminus = lld[j] + p;
#line 200 "dlaneg.f"
	    if (dminus < 0.) {
#line 200 "dlaneg.f"
		++neg2;
#line 200 "dlaneg.f"
	    }
#line 201 "dlaneg.f"
	    tmp = p / dminus;
#line 202 "dlaneg.f"
	    p = tmp * d__[j] - *sigma;
#line 203 "dlaneg.f"
/* L23: */
#line 203 "dlaneg.f"
	}
#line 204 "dlaneg.f"
	sawnan = disnan_(&p);
/*     As above, run a slower version that substitutes 1 for Inf/Inf. */

#line 207 "dlaneg.f"
	if (sawnan) {
#line 208 "dlaneg.f"
	    neg2 = 0;
#line 209 "dlaneg.f"
	    p = bsav;
/* Computing MAX */
#line 210 "dlaneg.f"
	    i__3 = bj - 127;
#line 210 "dlaneg.f"
	    i__2 = max(i__3,*r__);
#line 210 "dlaneg.f"
	    for (j = bj; j >= i__2; --j) {
#line 211 "dlaneg.f"
		dminus = lld[j] + p;
#line 212 "dlaneg.f"
		if (dminus < 0.) {
#line 212 "dlaneg.f"
		    ++neg2;
#line 212 "dlaneg.f"
		}
#line 213 "dlaneg.f"
		tmp = p / dminus;
#line 214 "dlaneg.f"
		if (disnan_(&tmp)) {
#line 214 "dlaneg.f"
		    tmp = 1.;
#line 214 "dlaneg.f"
		}
#line 215 "dlaneg.f"
		p = tmp * d__[j] - *sigma;
#line 216 "dlaneg.f"
/* L24: */
#line 216 "dlaneg.f"
	    }
#line 217 "dlaneg.f"
	}
#line 218 "dlaneg.f"
	negcnt += neg2;
#line 219 "dlaneg.f"
/* L230: */
#line 219 "dlaneg.f"
    }

/*     III) Twist index */
/*       T was shifted by SIGMA initially. */
#line 223 "dlaneg.f"
    gamma = t + *sigma + p;
#line 224 "dlaneg.f"
    if (gamma < 0.) {
#line 224 "dlaneg.f"
	++negcnt;
#line 224 "dlaneg.f"
    }
#line 226 "dlaneg.f"
    ret_val = negcnt;
#line 227 "dlaneg.f"
    return ret_val;
} /* dlaneg_ */

