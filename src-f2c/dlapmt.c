#line 1 "dlapmt.f"
/* dlapmt.f -- translated by f2c (version 20100827).
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

#line 1 "dlapmt.f"
/* > \brief \b DLAPMT performs a forward or backward permutation of the columns of a matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLAPMT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlapmt.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlapmt.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlapmt.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLAPMT( FORWRD, M, N, X, LDX, K ) */

/*       .. Scalar Arguments .. */
/*       LOGICAL            FORWRD */
/*       INTEGER            LDX, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            K( * ) */
/*       DOUBLE PRECISION   X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLAPMT rearranges the columns of the M by N matrix X as specified */
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
/* >          X is DOUBLE PRECISION array, dimension (LDX,N) */
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

/* > \ingroup doubleOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int dlapmt_(logical *forwrd, integer *m, integer *n, 
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

#line 129 "dlapmt.f"
    /* Parameter adjustments */
#line 129 "dlapmt.f"
    x_dim1 = *ldx;
#line 129 "dlapmt.f"
    x_offset = 1 + x_dim1;
#line 129 "dlapmt.f"
    x -= x_offset;
#line 129 "dlapmt.f"
    --k;
#line 129 "dlapmt.f"

#line 129 "dlapmt.f"
    /* Function Body */
#line 129 "dlapmt.f"
    if (*n <= 1) {
#line 129 "dlapmt.f"
	return 0;
#line 129 "dlapmt.f"
    }

#line 132 "dlapmt.f"
    i__1 = *n;
#line 132 "dlapmt.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 133 "dlapmt.f"
	k[i__] = -k[i__];
#line 134 "dlapmt.f"
/* L10: */
#line 134 "dlapmt.f"
    }

#line 136 "dlapmt.f"
    if (*forwrd) {

/*        Forward permutation */

#line 140 "dlapmt.f"
	i__1 = *n;
#line 140 "dlapmt.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

#line 142 "dlapmt.f"
	    if (k[i__] > 0) {
#line 142 "dlapmt.f"
		goto L40;
#line 142 "dlapmt.f"
	    }

#line 145 "dlapmt.f"
	    j = i__;
#line 146 "dlapmt.f"
	    k[j] = -k[j];
#line 147 "dlapmt.f"
	    in = k[j];

#line 149 "dlapmt.f"
L20:
#line 150 "dlapmt.f"
	    if (k[in] > 0) {
#line 150 "dlapmt.f"
		goto L40;
#line 150 "dlapmt.f"
	    }

#line 153 "dlapmt.f"
	    i__2 = *m;
#line 153 "dlapmt.f"
	    for (ii = 1; ii <= i__2; ++ii) {
#line 154 "dlapmt.f"
		temp = x[ii + j * x_dim1];
#line 155 "dlapmt.f"
		x[ii + j * x_dim1] = x[ii + in * x_dim1];
#line 156 "dlapmt.f"
		x[ii + in * x_dim1] = temp;
#line 157 "dlapmt.f"
/* L30: */
#line 157 "dlapmt.f"
	    }

#line 159 "dlapmt.f"
	    k[in] = -k[in];
#line 160 "dlapmt.f"
	    j = in;
#line 161 "dlapmt.f"
	    in = k[in];
#line 162 "dlapmt.f"
	    goto L20;

#line 164 "dlapmt.f"
L40:

#line 166 "dlapmt.f"
/* L50: */
#line 166 "dlapmt.f"
	    ;
#line 166 "dlapmt.f"
	}

#line 168 "dlapmt.f"
    } else {

/*        Backward permutation */

#line 172 "dlapmt.f"
	i__1 = *n;
#line 172 "dlapmt.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

#line 174 "dlapmt.f"
	    if (k[i__] > 0) {
#line 174 "dlapmt.f"
		goto L80;
#line 174 "dlapmt.f"
	    }

#line 177 "dlapmt.f"
	    k[i__] = -k[i__];
#line 178 "dlapmt.f"
	    j = k[i__];
#line 179 "dlapmt.f"
L60:
#line 180 "dlapmt.f"
	    if (j == i__) {
#line 180 "dlapmt.f"
		goto L80;
#line 180 "dlapmt.f"
	    }

#line 183 "dlapmt.f"
	    i__2 = *m;
#line 183 "dlapmt.f"
	    for (ii = 1; ii <= i__2; ++ii) {
#line 184 "dlapmt.f"
		temp = x[ii + i__ * x_dim1];
#line 185 "dlapmt.f"
		x[ii + i__ * x_dim1] = x[ii + j * x_dim1];
#line 186 "dlapmt.f"
		x[ii + j * x_dim1] = temp;
#line 187 "dlapmt.f"
/* L70: */
#line 187 "dlapmt.f"
	    }

#line 189 "dlapmt.f"
	    k[j] = -k[j];
#line 190 "dlapmt.f"
	    j = k[j];
#line 191 "dlapmt.f"
	    goto L60;

#line 193 "dlapmt.f"
L80:

#line 195 "dlapmt.f"
/* L90: */
#line 195 "dlapmt.f"
	    ;
#line 195 "dlapmt.f"
	}

#line 197 "dlapmt.f"
    }

#line 199 "dlapmt.f"
    return 0;

/*     End of DLAPMT */

} /* dlapmt_ */

