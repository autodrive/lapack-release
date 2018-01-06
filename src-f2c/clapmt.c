#line 1 "clapmt.f"
/* clapmt.f -- translated by f2c (version 20100827).
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

#line 1 "clapmt.f"
/* > \brief \b CLAPMT performs a forward or backward permutation of the columns of a matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLAPMT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clapmt.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clapmt.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clapmt.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLAPMT( FORWRD, M, N, X, LDX, K ) */

/*       .. Scalar Arguments .. */
/*       LOGICAL            FORWRD */
/*       INTEGER            LDX, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            K( * ) */
/*       COMPLEX            X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLAPMT rearranges the columns of the M by N matrix X as specified */
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
/* >          X is COMPLEX array, dimension (LDX,N) */
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

/* > \ingroup complexOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int clapmt_(logical *forwrd, integer *m, integer *n, 
	doublecomplex *x, integer *ldx, integer *k)
{
    /* System generated locals */
    integer x_dim1, x_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, j, ii, in;
    static doublecomplex temp;


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

#line 129 "clapmt.f"
    /* Parameter adjustments */
#line 129 "clapmt.f"
    x_dim1 = *ldx;
#line 129 "clapmt.f"
    x_offset = 1 + x_dim1;
#line 129 "clapmt.f"
    x -= x_offset;
#line 129 "clapmt.f"
    --k;
#line 129 "clapmt.f"

#line 129 "clapmt.f"
    /* Function Body */
#line 129 "clapmt.f"
    if (*n <= 1) {
#line 129 "clapmt.f"
	return 0;
#line 129 "clapmt.f"
    }

#line 132 "clapmt.f"
    i__1 = *n;
#line 132 "clapmt.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 133 "clapmt.f"
	k[i__] = -k[i__];
#line 134 "clapmt.f"
/* L10: */
#line 134 "clapmt.f"
    }

#line 136 "clapmt.f"
    if (*forwrd) {

/*        Forward permutation */

#line 140 "clapmt.f"
	i__1 = *n;
#line 140 "clapmt.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

#line 142 "clapmt.f"
	    if (k[i__] > 0) {
#line 142 "clapmt.f"
		goto L40;
#line 142 "clapmt.f"
	    }

#line 145 "clapmt.f"
	    j = i__;
#line 146 "clapmt.f"
	    k[j] = -k[j];
#line 147 "clapmt.f"
	    in = k[j];

#line 149 "clapmt.f"
L20:
#line 150 "clapmt.f"
	    if (k[in] > 0) {
#line 150 "clapmt.f"
		goto L40;
#line 150 "clapmt.f"
	    }

#line 153 "clapmt.f"
	    i__2 = *m;
#line 153 "clapmt.f"
	    for (ii = 1; ii <= i__2; ++ii) {
#line 154 "clapmt.f"
		i__3 = ii + j * x_dim1;
#line 154 "clapmt.f"
		temp.r = x[i__3].r, temp.i = x[i__3].i;
#line 155 "clapmt.f"
		i__3 = ii + j * x_dim1;
#line 155 "clapmt.f"
		i__4 = ii + in * x_dim1;
#line 155 "clapmt.f"
		x[i__3].r = x[i__4].r, x[i__3].i = x[i__4].i;
#line 156 "clapmt.f"
		i__3 = ii + in * x_dim1;
#line 156 "clapmt.f"
		x[i__3].r = temp.r, x[i__3].i = temp.i;
#line 157 "clapmt.f"
/* L30: */
#line 157 "clapmt.f"
	    }

#line 159 "clapmt.f"
	    k[in] = -k[in];
#line 160 "clapmt.f"
	    j = in;
#line 161 "clapmt.f"
	    in = k[in];
#line 162 "clapmt.f"
	    goto L20;

#line 164 "clapmt.f"
L40:

#line 166 "clapmt.f"
/* L60: */
#line 166 "clapmt.f"
	    ;
#line 166 "clapmt.f"
	}

#line 168 "clapmt.f"
    } else {

/*        Backward permutation */

#line 172 "clapmt.f"
	i__1 = *n;
#line 172 "clapmt.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

#line 174 "clapmt.f"
	    if (k[i__] > 0) {
#line 174 "clapmt.f"
		goto L100;
#line 174 "clapmt.f"
	    }

#line 177 "clapmt.f"
	    k[i__] = -k[i__];
#line 178 "clapmt.f"
	    j = k[i__];
#line 179 "clapmt.f"
L80:
#line 180 "clapmt.f"
	    if (j == i__) {
#line 180 "clapmt.f"
		goto L100;
#line 180 "clapmt.f"
	    }

#line 183 "clapmt.f"
	    i__2 = *m;
#line 183 "clapmt.f"
	    for (ii = 1; ii <= i__2; ++ii) {
#line 184 "clapmt.f"
		i__3 = ii + i__ * x_dim1;
#line 184 "clapmt.f"
		temp.r = x[i__3].r, temp.i = x[i__3].i;
#line 185 "clapmt.f"
		i__3 = ii + i__ * x_dim1;
#line 185 "clapmt.f"
		i__4 = ii + j * x_dim1;
#line 185 "clapmt.f"
		x[i__3].r = x[i__4].r, x[i__3].i = x[i__4].i;
#line 186 "clapmt.f"
		i__3 = ii + j * x_dim1;
#line 186 "clapmt.f"
		x[i__3].r = temp.r, x[i__3].i = temp.i;
#line 187 "clapmt.f"
/* L90: */
#line 187 "clapmt.f"
	    }

#line 189 "clapmt.f"
	    k[j] = -k[j];
#line 190 "clapmt.f"
	    j = k[j];
#line 191 "clapmt.f"
	    goto L80;

#line 193 "clapmt.f"
L100:
#line 195 "clapmt.f"
/* L110: */
#line 195 "clapmt.f"
	    ;
#line 195 "clapmt.f"
	}

#line 197 "clapmt.f"
    }

#line 199 "clapmt.f"
    return 0;

/*     End of CLAPMT */

} /* clapmt_ */

