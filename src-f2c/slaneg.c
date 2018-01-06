#line 1 "slaneg.f"
/* slaneg.f -- translated by f2c (version 20100827).
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

#line 1 "slaneg.f"
/* > \brief \b SLANEG computes the Sturm count. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLANEG + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaneg.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaneg.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaneg.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       INTEGER FUNCTION SLANEG( N, D, LLD, SIGMA, PIVMIN, R ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            N, R */
/*       REAL               PIVMIN, SIGMA */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               D( * ), LLD( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLANEG computes the Sturm count, the number of negative pivots */
/* > encountered while factoring tridiagonal T - sigma I = L D L^T. */
/* > This implementation works directly on the factors without forming */
/* > the tridiagonal matrix T.  The Sturm count is also the number of */
/* > eigenvalues of T less than sigma. */
/* > */
/* > This routine is called from SLARRB. */
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
/* >          D is REAL array, dimension (N) */
/* >          The N diagonal elements of the diagonal matrix D. */
/* > \endverbatim */
/* > */
/* > \param[in] LLD */
/* > \verbatim */
/* >          LLD is REAL array, dimension (N-1) */
/* >          The (N-1) elements L(i)*L(i)*D(i). */
/* > \endverbatim */
/* > */
/* > \param[in] SIGMA */
/* > \verbatim */
/* >          SIGMA is REAL */
/* >          Shift amount in T - sigma I = L D L^T. */
/* > \endverbatim */
/* > */
/* > \param[in] PIVMIN */
/* > \verbatim */
/* >          PIVMIN is REAL */
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

/* > \date September 2012 */

/* > \ingroup auxOTHERauxiliary */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Osni Marques, LBNL/NERSC, USA \n */
/* >     Christof Voemel, University of California, Berkeley, USA \n */
/* >     Jason Riedy, University of California, Berkeley, USA \n */
/* > */
/*  ===================================================================== */
integer slaneg_(integer *n, doublereal *d__, doublereal *lld, doublereal *
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
    static integer negcnt;
    static logical sawnan;
    extern logical sisnan_(doublereal *);
    static doublereal dminus;


/*  -- LAPACK auxiliary routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

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
#line 161 "slaneg.f"
    /* Parameter adjustments */
#line 161 "slaneg.f"
    --lld;
#line 161 "slaneg.f"
    --d__;
#line 161 "slaneg.f"

#line 161 "slaneg.f"
    /* Function Body */
#line 161 "slaneg.f"
    negcnt = 0;
/*     I) upper part: L D L^T - SIGMA I = L+ D+ L+^T */
#line 164 "slaneg.f"
    t = -(*sigma);
#line 165 "slaneg.f"
    i__1 = *r__ - 1;
#line 165 "slaneg.f"
    for (bj = 1; bj <= i__1; bj += 128) {
#line 166 "slaneg.f"
	neg1 = 0;
#line 167 "slaneg.f"
	bsav = t;
/* Computing MIN */
#line 168 "slaneg.f"
	i__3 = bj + 127, i__4 = *r__ - 1;
#line 168 "slaneg.f"
	i__2 = min(i__3,i__4);
#line 168 "slaneg.f"
	for (j = bj; j <= i__2; ++j) {
#line 169 "slaneg.f"
	    dplus = d__[j] + t;
#line 170 "slaneg.f"
	    if (dplus < 0.) {
#line 170 "slaneg.f"
		++neg1;
#line 170 "slaneg.f"
	    }
#line 171 "slaneg.f"
	    tmp = t / dplus;
#line 172 "slaneg.f"
	    t = tmp * lld[j] - *sigma;
#line 173 "slaneg.f"
/* L21: */
#line 173 "slaneg.f"
	}
#line 174 "slaneg.f"
	sawnan = sisnan_(&t);
/*     Run a slower version of the above loop if a NaN is detected. */
/*     A NaN should occur only with a zero pivot after an infinite */
/*     pivot.  In that case, substituting 1 for T/DPLUS is the */
/*     correct limit. */
#line 179 "slaneg.f"
	if (sawnan) {
#line 180 "slaneg.f"
	    neg1 = 0;
#line 181 "slaneg.f"
	    t = bsav;
/* Computing MIN */
#line 182 "slaneg.f"
	    i__3 = bj + 127, i__4 = *r__ - 1;
#line 182 "slaneg.f"
	    i__2 = min(i__3,i__4);
#line 182 "slaneg.f"
	    for (j = bj; j <= i__2; ++j) {
#line 183 "slaneg.f"
		dplus = d__[j] + t;
#line 184 "slaneg.f"
		if (dplus < 0.) {
#line 184 "slaneg.f"
		    ++neg1;
#line 184 "slaneg.f"
		}
#line 185 "slaneg.f"
		tmp = t / dplus;
#line 186 "slaneg.f"
		if (sisnan_(&tmp)) {
#line 186 "slaneg.f"
		    tmp = 1.;
#line 186 "slaneg.f"
		}
#line 187 "slaneg.f"
		t = tmp * lld[j] - *sigma;
#line 188 "slaneg.f"
/* L22: */
#line 188 "slaneg.f"
	    }
#line 189 "slaneg.f"
	}
#line 190 "slaneg.f"
	negcnt += neg1;
#line 191 "slaneg.f"
/* L210: */
#line 191 "slaneg.f"
    }

/*     II) lower part: L D L^T - SIGMA I = U- D- U-^T */
#line 194 "slaneg.f"
    p = d__[*n] - *sigma;
#line 195 "slaneg.f"
    i__1 = *r__;
#line 195 "slaneg.f"
    for (bj = *n - 1; bj >= i__1; bj += -128) {
#line 196 "slaneg.f"
	neg2 = 0;
#line 197 "slaneg.f"
	bsav = p;
/* Computing MAX */
#line 198 "slaneg.f"
	i__3 = bj - 127;
#line 198 "slaneg.f"
	i__2 = max(i__3,*r__);
#line 198 "slaneg.f"
	for (j = bj; j >= i__2; --j) {
#line 199 "slaneg.f"
	    dminus = lld[j] + p;
#line 200 "slaneg.f"
	    if (dminus < 0.) {
#line 200 "slaneg.f"
		++neg2;
#line 200 "slaneg.f"
	    }
#line 201 "slaneg.f"
	    tmp = p / dminus;
#line 202 "slaneg.f"
	    p = tmp * d__[j] - *sigma;
#line 203 "slaneg.f"
/* L23: */
#line 203 "slaneg.f"
	}
#line 204 "slaneg.f"
	sawnan = sisnan_(&p);
/*     As above, run a slower version that substitutes 1 for Inf/Inf. */

#line 207 "slaneg.f"
	if (sawnan) {
#line 208 "slaneg.f"
	    neg2 = 0;
#line 209 "slaneg.f"
	    p = bsav;
/* Computing MAX */
#line 210 "slaneg.f"
	    i__3 = bj - 127;
#line 210 "slaneg.f"
	    i__2 = max(i__3,*r__);
#line 210 "slaneg.f"
	    for (j = bj; j >= i__2; --j) {
#line 211 "slaneg.f"
		dminus = lld[j] + p;
#line 212 "slaneg.f"
		if (dminus < 0.) {
#line 212 "slaneg.f"
		    ++neg2;
#line 212 "slaneg.f"
		}
#line 213 "slaneg.f"
		tmp = p / dminus;
#line 214 "slaneg.f"
		if (sisnan_(&tmp)) {
#line 214 "slaneg.f"
		    tmp = 1.;
#line 214 "slaneg.f"
		}
#line 215 "slaneg.f"
		p = tmp * d__[j] - *sigma;
#line 216 "slaneg.f"
/* L24: */
#line 216 "slaneg.f"
	    }
#line 217 "slaneg.f"
	}
#line 218 "slaneg.f"
	negcnt += neg2;
#line 219 "slaneg.f"
/* L230: */
#line 219 "slaneg.f"
    }

/*     III) Twist index */
/*       T was shifted by SIGMA initially. */
#line 223 "slaneg.f"
    gamma = t + *sigma + p;
#line 224 "slaneg.f"
    if (gamma < 0.) {
#line 224 "slaneg.f"
	++negcnt;
#line 224 "slaneg.f"
    }
#line 226 "slaneg.f"
    ret_val = negcnt;
#line 227 "slaneg.f"
    return ret_val;
} /* slaneg_ */

