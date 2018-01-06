#line 1 "slapmt.f"
/* slapmt.f -- translated by f2c (version 20100827).
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

#line 1 "slapmt.f"
/* > \brief \b SLAPMT performs a forward or backward permutation of the columns of a matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLAPMT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slapmt.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slapmt.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slapmt.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLAPMT( FORWRD, M, N, X, LDX, K ) */

/*       .. Scalar Arguments .. */
/*       LOGICAL            FORWRD */
/*       INTEGER            LDX, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            K( * ) */
/*       REAL               X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLAPMT rearranges the columns of the M by N matrix X as specified */
/* > by the permutation K(1),K(2),...,K(N) of the integers 1,...,N. */
/* > If FORWRD = .TRUE.,  forward permutation: */
/* > */
/* >      X(*,K(J)) is moved X(*,J) for J = 1,2,...,N. */
/* > */
/* > If FORWRD = .FALSE., backward permutation: */
/* > */
/* >      X(*,J) is moved to X(*,K(J)) for J = 1,2,...,N. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] FORWRD */
/* > \verbatim */
/* >          FORWRD is LOGICAL */
/* >          = .TRUE., forward permutation */
/* >          = .FALSE., backward permutation */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix X. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix X. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* >          X is REAL array, dimension (LDX,N) */
/* >          On entry, the M by N matrix X. */
/* >          On exit, X contains the permuted matrix X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX */
/* > \verbatim */
/* >          LDX is INTEGER */
/* >          The leading dimension of the array X, LDX >= MAX(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in,out] K */
/* > \verbatim */
/* >          K is INTEGER array, dimension (N) */
/* >          On entry, K contains the permutation vector. K is used as */
/* >          internal workspace, but reset to its original value on */
/* >          output. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup realOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int slapmt_(logical *forwrd, integer *m, integer *n, 
	doublereal *x, integer *ldx, integer *k)
{
    /* System generated locals */
    integer x_dim1, x_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, ii, in;
    static doublereal temp;


/*  -- LAPACK auxiliary routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */

#line 129 "slapmt.f"
    /* Parameter adjustments */
#line 129 "slapmt.f"
    x_dim1 = *ldx;
#line 129 "slapmt.f"
    x_offset = 1 + x_dim1;
#line 129 "slapmt.f"
    x -= x_offset;
#line 129 "slapmt.f"
    --k;
#line 129 "slapmt.f"

#line 129 "slapmt.f"
    /* Function Body */
#line 129 "slapmt.f"
    if (*n <= 1) {
#line 129 "slapmt.f"
	return 0;
#line 129 "slapmt.f"
    }

#line 132 "slapmt.f"
    i__1 = *n;
#line 132 "slapmt.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 133 "slapmt.f"
	k[i__] = -k[i__];
#line 134 "slapmt.f"
/* L10: */
#line 134 "slapmt.f"
    }

#line 136 "slapmt.f"
    if (*forwrd) {

/*        Forward permutation */

#line 140 "slapmt.f"
	i__1 = *n;
#line 140 "slapmt.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

#line 142 "slapmt.f"
	    if (k[i__] > 0) {
#line 142 "slapmt.f"
		goto L40;
#line 142 "slapmt.f"
	    }

#line 145 "slapmt.f"
	    j = i__;
#line 146 "slapmt.f"
	    k[j] = -k[j];
#line 147 "slapmt.f"
	    in = k[j];

#line 149 "slapmt.f"
L20:
#line 150 "slapmt.f"
	    if (k[in] > 0) {
#line 150 "slapmt.f"
		goto L40;
#line 150 "slapmt.f"
	    }

#line 153 "slapmt.f"
	    i__2 = *m;
#line 153 "slapmt.f"
	    for (ii = 1; ii <= i__2; ++ii) {
#line 154 "slapmt.f"
		temp = x[ii + j * x_dim1];
#line 155 "slapmt.f"
		x[ii + j * x_dim1] = x[ii + in * x_dim1];
#line 156 "slapmt.f"
		x[ii + in * x_dim1] = temp;
#line 157 "slapmt.f"
/* L30: */
#line 157 "slapmt.f"
	    }

#line 159 "slapmt.f"
	    k[in] = -k[in];
#line 160 "slapmt.f"
	    j = in;
#line 161 "slapmt.f"
	    in = k[in];
#line 162 "slapmt.f"
	    goto L20;

#line 164 "slapmt.f"
L40:

#line 166 "slapmt.f"
/* L60: */
#line 166 "slapmt.f"
	    ;
#line 166 "slapmt.f"
	}

#line 168 "slapmt.f"
    } else {

/*        Backward permutation */

#line 172 "slapmt.f"
	i__1 = *n;
#line 172 "slapmt.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

#line 174 "slapmt.f"
	    if (k[i__] > 0) {
#line 174 "slapmt.f"
		goto L100;
#line 174 "slapmt.f"
	    }

#line 177 "slapmt.f"
	    k[i__] = -k[i__];
#line 178 "slapmt.f"
	    j = k[i__];
#line 179 "slapmt.f"
L80:
#line 180 "slapmt.f"
	    if (j == i__) {
#line 180 "slapmt.f"
		goto L100;
#line 180 "slapmt.f"
	    }

#line 183 "slapmt.f"
	    i__2 = *m;
#line 183 "slapmt.f"
	    for (ii = 1; ii <= i__2; ++ii) {
#line 184 "slapmt.f"
		temp = x[ii + i__ * x_dim1];
#line 185 "slapmt.f"
		x[ii + i__ * x_dim1] = x[ii + j * x_dim1];
#line 186 "slapmt.f"
		x[ii + j * x_dim1] = temp;
#line 187 "slapmt.f"
/* L90: */
#line 187 "slapmt.f"
	    }

#line 189 "slapmt.f"
	    k[j] = -k[j];
#line 190 "slapmt.f"
	    j = k[j];
#line 191 "slapmt.f"
	    goto L80;

#line 193 "slapmt.f"
L100:
#line 195 "slapmt.f"
/* L110: */
#line 195 "slapmt.f"
	    ;
#line 195 "slapmt.f"
	}

#line 197 "slapmt.f"
    }

#line 199 "slapmt.f"
    return 0;

/*     End of SLAPMT */

} /* slapmt_ */

