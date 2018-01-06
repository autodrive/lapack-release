#line 1 "zlapmt.f"
/* zlapmt.f -- translated by f2c (version 20100827).
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

#line 1 "zlapmt.f"
/* > \brief \b ZLAPMT performs a forward or backward permutation of the columns of a matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLAPMT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlapmt.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlapmt.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlapmt.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLAPMT( FORWRD, M, N, X, LDX, K ) */

/*       .. Scalar Arguments .. */
/*       LOGICAL            FORWRD */
/*       INTEGER            LDX, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            K( * ) */
/*       COMPLEX*16         X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLAPMT rearranges the columns of the M by N matrix X as specified */
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
/* >          X is COMPLEX*16 array, dimension (LDX,N) */
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

/* > \ingroup complex16OTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int zlapmt_(logical *forwrd, integer *m, integer *n, 
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

#line 129 "zlapmt.f"
    /* Parameter adjustments */
#line 129 "zlapmt.f"
    x_dim1 = *ldx;
#line 129 "zlapmt.f"
    x_offset = 1 + x_dim1;
#line 129 "zlapmt.f"
    x -= x_offset;
#line 129 "zlapmt.f"
    --k;
#line 129 "zlapmt.f"

#line 129 "zlapmt.f"
    /* Function Body */
#line 129 "zlapmt.f"
    if (*n <= 1) {
#line 129 "zlapmt.f"
	return 0;
#line 129 "zlapmt.f"
    }

#line 132 "zlapmt.f"
    i__1 = *n;
#line 132 "zlapmt.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 133 "zlapmt.f"
	k[i__] = -k[i__];
#line 134 "zlapmt.f"
/* L10: */
#line 134 "zlapmt.f"
    }

#line 136 "zlapmt.f"
    if (*forwrd) {

/*        Forward permutation */

#line 140 "zlapmt.f"
	i__1 = *n;
#line 140 "zlapmt.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

#line 142 "zlapmt.f"
	    if (k[i__] > 0) {
#line 142 "zlapmt.f"
		goto L40;
#line 142 "zlapmt.f"
	    }

#line 145 "zlapmt.f"
	    j = i__;
#line 146 "zlapmt.f"
	    k[j] = -k[j];
#line 147 "zlapmt.f"
	    in = k[j];

#line 149 "zlapmt.f"
L20:
#line 150 "zlapmt.f"
	    if (k[in] > 0) {
#line 150 "zlapmt.f"
		goto L40;
#line 150 "zlapmt.f"
	    }

#line 153 "zlapmt.f"
	    i__2 = *m;
#line 153 "zlapmt.f"
	    for (ii = 1; ii <= i__2; ++ii) {
#line 154 "zlapmt.f"
		i__3 = ii + j * x_dim1;
#line 154 "zlapmt.f"
		temp.r = x[i__3].r, temp.i = x[i__3].i;
#line 155 "zlapmt.f"
		i__3 = ii + j * x_dim1;
#line 155 "zlapmt.f"
		i__4 = ii + in * x_dim1;
#line 155 "zlapmt.f"
		x[i__3].r = x[i__4].r, x[i__3].i = x[i__4].i;
#line 156 "zlapmt.f"
		i__3 = ii + in * x_dim1;
#line 156 "zlapmt.f"
		x[i__3].r = temp.r, x[i__3].i = temp.i;
#line 157 "zlapmt.f"
/* L30: */
#line 157 "zlapmt.f"
	    }

#line 159 "zlapmt.f"
	    k[in] = -k[in];
#line 160 "zlapmt.f"
	    j = in;
#line 161 "zlapmt.f"
	    in = k[in];
#line 162 "zlapmt.f"
	    goto L20;

#line 164 "zlapmt.f"
L40:

#line 166 "zlapmt.f"
/* L50: */
#line 166 "zlapmt.f"
	    ;
#line 166 "zlapmt.f"
	}

#line 168 "zlapmt.f"
    } else {

/*        Backward permutation */

#line 172 "zlapmt.f"
	i__1 = *n;
#line 172 "zlapmt.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

#line 174 "zlapmt.f"
	    if (k[i__] > 0) {
#line 174 "zlapmt.f"
		goto L80;
#line 174 "zlapmt.f"
	    }

#line 177 "zlapmt.f"
	    k[i__] = -k[i__];
#line 178 "zlapmt.f"
	    j = k[i__];
#line 179 "zlapmt.f"
L60:
#line 180 "zlapmt.f"
	    if (j == i__) {
#line 180 "zlapmt.f"
		goto L80;
#line 180 "zlapmt.f"
	    }

#line 183 "zlapmt.f"
	    i__2 = *m;
#line 183 "zlapmt.f"
	    for (ii = 1; ii <= i__2; ++ii) {
#line 184 "zlapmt.f"
		i__3 = ii + i__ * x_dim1;
#line 184 "zlapmt.f"
		temp.r = x[i__3].r, temp.i = x[i__3].i;
#line 185 "zlapmt.f"
		i__3 = ii + i__ * x_dim1;
#line 185 "zlapmt.f"
		i__4 = ii + j * x_dim1;
#line 185 "zlapmt.f"
		x[i__3].r = x[i__4].r, x[i__3].i = x[i__4].i;
#line 186 "zlapmt.f"
		i__3 = ii + j * x_dim1;
#line 186 "zlapmt.f"
		x[i__3].r = temp.r, x[i__3].i = temp.i;
#line 187 "zlapmt.f"
/* L70: */
#line 187 "zlapmt.f"
	    }

#line 189 "zlapmt.f"
	    k[j] = -k[j];
#line 190 "zlapmt.f"
	    j = k[j];
#line 191 "zlapmt.f"
	    goto L60;

#line 193 "zlapmt.f"
L80:

#line 195 "zlapmt.f"
/* L90: */
#line 195 "zlapmt.f"
	    ;
#line 195 "zlapmt.f"
	}

#line 197 "zlapmt.f"
    }

#line 199 "zlapmt.f"
    return 0;

/*     End of ZLAPMT */

} /* zlapmt_ */

