#line 1 "dlarra.f"
/* dlarra.f -- translated by f2c (version 20100827).
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

#line 1 "dlarra.f"
/* > \brief \b DLARRA computes the splitting points with the specified threshold. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLARRA + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarra.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarra.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarra.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLARRA( N, D, E, E2, SPLTOL, TNRM, */
/*                           NSPLIT, ISPLIT, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, N, NSPLIT */
/*       DOUBLE PRECISION    SPLTOL, TNRM */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            ISPLIT( * ) */
/*       DOUBLE PRECISION   D( * ), E( * ), E2( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > Compute the splitting points with threshold SPLTOL. */
/* > DLARRA sets any "small" off-diagonal elements to zero. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix. N > 0. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          On entry, the N diagonal elements of the tridiagonal */
/* >          matrix T. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array, dimension (N) */
/* >          On entry, the first (N-1) entries contain the subdiagonal */
/* >          elements of the tridiagonal matrix T; E(N) need not be set. */
/* >          On exit, the entries E( ISPLIT( I ) ), 1 <= I <= NSPLIT, */
/* >          are set to zero, the other entries of E are untouched. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E2 */
/* > \verbatim */
/* >          E2 is DOUBLE PRECISION array, dimension (N) */
/* >          On entry, the first (N-1) entries contain the SQUARES of the */
/* >          subdiagonal elements of the tridiagonal matrix T; */
/* >          E2(N) need not be set. */
/* >          On exit, the entries E2( ISPLIT( I ) ), */
/* >          1 <= I <= NSPLIT, have been set to zero */
/* > \endverbatim */
/* > */
/* > \param[in] SPLTOL */
/* > \verbatim */
/* >          SPLTOL is DOUBLE PRECISION */
/* >          The threshold for splitting. Two criteria can be used: */
/* >          SPLTOL<0 : criterion based on absolute off-diagonal value */
/* >          SPLTOL>0 : criterion that preserves relative accuracy */
/* > \endverbatim */
/* > */
/* > \param[in] TNRM */
/* > \verbatim */
/* >          TNRM is DOUBLE PRECISION */
/* >          The norm of the matrix. */
/* > \endverbatim */
/* > */
/* > \param[out] NSPLIT */
/* > \verbatim */
/* >          NSPLIT is INTEGER */
/* >          The number of blocks T splits into. 1 <= NSPLIT <= N. */
/* > \endverbatim */
/* > */
/* > \param[out] ISPLIT */
/* > \verbatim */
/* >          ISPLIT is INTEGER array, dimension (N) */
/* >          The splitting points, at which T breaks up into blocks. */
/* >          The first block consists of rows/columns 1 to ISPLIT(1), */
/* >          the second of rows/columns ISPLIT(1)+1 through ISPLIT(2), */
/* >          etc., and the NSPLIT-th consists of rows/columns */
/* >          ISPLIT(NSPLIT-1)+1 through ISPLIT(NSPLIT)=N. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
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
/* > Beresford Parlett, University of California, Berkeley, USA \n */
/* > Jim Demmel, University of California, Berkeley, USA \n */
/* > Inderjit Dhillon, University of Texas, Austin, USA \n */
/* > Osni Marques, LBNL/NERSC, USA \n */
/* > Christof Voemel, University of California, Berkeley, USA */

/*  ===================================================================== */
/* Subroutine */ int dlarra_(integer *n, doublereal *d__, doublereal *e, 
	doublereal *e2, doublereal *spltol, doublereal *tnrm, integer *nsplit,
	 integer *isplit, integer *info)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal tmp1, eabs;


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
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 169 "dlarra.f"
    /* Parameter adjustments */
#line 169 "dlarra.f"
    --isplit;
#line 169 "dlarra.f"
    --e2;
#line 169 "dlarra.f"
    --e;
#line 169 "dlarra.f"
    --d__;
#line 169 "dlarra.f"

#line 169 "dlarra.f"
    /* Function Body */
#line 169 "dlarra.f"
    *info = 0;
/*     Compute splitting points */
#line 172 "dlarra.f"
    *nsplit = 1;
#line 173 "dlarra.f"
    if (*spltol < 0.) {
/*        Criterion based on absolute off-diagonal value */
#line 175 "dlarra.f"
	tmp1 = abs(*spltol) * *tnrm;
#line 176 "dlarra.f"
	i__1 = *n - 1;
#line 176 "dlarra.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 177 "dlarra.f"
	    eabs = (d__1 = e[i__], abs(d__1));
#line 178 "dlarra.f"
	    if (eabs <= tmp1) {
#line 179 "dlarra.f"
		e[i__] = 0.;
#line 180 "dlarra.f"
		e2[i__] = 0.;
#line 181 "dlarra.f"
		isplit[*nsplit] = i__;
#line 182 "dlarra.f"
		++(*nsplit);
#line 183 "dlarra.f"
	    }
#line 184 "dlarra.f"
/* L9: */
#line 184 "dlarra.f"
	}
#line 185 "dlarra.f"
    } else {
/*        Criterion that guarantees relative accuracy */
#line 187 "dlarra.f"
	i__1 = *n - 1;
#line 187 "dlarra.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 188 "dlarra.f"
	    eabs = (d__1 = e[i__], abs(d__1));
#line 189 "dlarra.f"
	    if (eabs <= *spltol * sqrt((d__1 = d__[i__], abs(d__1))) * sqrt((
		    d__2 = d__[i__ + 1], abs(d__2)))) {
#line 191 "dlarra.f"
		e[i__] = 0.;
#line 192 "dlarra.f"
		e2[i__] = 0.;
#line 193 "dlarra.f"
		isplit[*nsplit] = i__;
#line 194 "dlarra.f"
		++(*nsplit);
#line 195 "dlarra.f"
	    }
#line 196 "dlarra.f"
/* L10: */
#line 196 "dlarra.f"
	}
#line 197 "dlarra.f"
    }
#line 198 "dlarra.f"
    isplit[*nsplit] = *n;
#line 200 "dlarra.f"
    return 0;

/*     End of DLARRA */

} /* dlarra_ */

