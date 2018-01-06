#line 1 "clapmr.f"
/* clapmr.f -- translated by f2c (version 20100827).
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

#line 1 "clapmr.f"
/* > \brief \b CLAPMR rearranges rows of a matrix as specified by a permutation vector. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLAPMR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clapmr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clapmr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clapmr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLAPMR( FORWRD, M, N, X, LDX, K ) */

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
/* > CLAPMR rearranges the rows of the M by N matrix X as specified */
/* > by the permutation K(1),K(2),...,K(M) of the integers 1,...,M. */
/* > If FORWRD = .TRUE.,  forward permutation: */
/* > */
/* >      X(K(I),*) is moved X(I,*) for I = 1,2,...,M. */
/* > */
/* > If FORWRD = .FALSE., backward permutation: */
/* > */
/* >      X(I,*) is moved to X(K(I),*) for I = 1,2,...,M. */
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
/* >          K is INTEGER array, dimension (M) */
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

/* > \date December 2016 */

/* > \ingroup complexOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int clapmr_(logical *forwrd, integer *m, integer *n, 
	doublecomplex *x, integer *ldx, integer *k)
{
    /* System generated locals */
    integer x_dim1, x_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, j, jj, in;
    static doublecomplex temp;


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */

#line 129 "clapmr.f"
    /* Parameter adjustments */
#line 129 "clapmr.f"
    x_dim1 = *ldx;
#line 129 "clapmr.f"
    x_offset = 1 + x_dim1;
#line 129 "clapmr.f"
    x -= x_offset;
#line 129 "clapmr.f"
    --k;
#line 129 "clapmr.f"

#line 129 "clapmr.f"
    /* Function Body */
#line 129 "clapmr.f"
    if (*m <= 1) {
#line 129 "clapmr.f"
	return 0;
#line 129 "clapmr.f"
    }

#line 132 "clapmr.f"
    i__1 = *m;
#line 132 "clapmr.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 133 "clapmr.f"
	k[i__] = -k[i__];
#line 134 "clapmr.f"
/* L10: */
#line 134 "clapmr.f"
    }

#line 136 "clapmr.f"
    if (*forwrd) {

/*        Forward permutation */

#line 140 "clapmr.f"
	i__1 = *m;
#line 140 "clapmr.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

#line 142 "clapmr.f"
	    if (k[i__] > 0) {
#line 142 "clapmr.f"
		goto L40;
#line 142 "clapmr.f"
	    }

#line 145 "clapmr.f"
	    j = i__;
#line 146 "clapmr.f"
	    k[j] = -k[j];
#line 147 "clapmr.f"
	    in = k[j];

#line 149 "clapmr.f"
L20:
#line 150 "clapmr.f"
	    if (k[in] > 0) {
#line 150 "clapmr.f"
		goto L40;
#line 150 "clapmr.f"
	    }

#line 153 "clapmr.f"
	    i__2 = *n;
#line 153 "clapmr.f"
	    for (jj = 1; jj <= i__2; ++jj) {
#line 154 "clapmr.f"
		i__3 = j + jj * x_dim1;
#line 154 "clapmr.f"
		temp.r = x[i__3].r, temp.i = x[i__3].i;
#line 155 "clapmr.f"
		i__3 = j + jj * x_dim1;
#line 155 "clapmr.f"
		i__4 = in + jj * x_dim1;
#line 155 "clapmr.f"
		x[i__3].r = x[i__4].r, x[i__3].i = x[i__4].i;
#line 156 "clapmr.f"
		i__3 = in + jj * x_dim1;
#line 156 "clapmr.f"
		x[i__3].r = temp.r, x[i__3].i = temp.i;
#line 157 "clapmr.f"
/* L30: */
#line 157 "clapmr.f"
	    }

#line 159 "clapmr.f"
	    k[in] = -k[in];
#line 160 "clapmr.f"
	    j = in;
#line 161 "clapmr.f"
	    in = k[in];
#line 162 "clapmr.f"
	    goto L20;

#line 164 "clapmr.f"
L40:

#line 166 "clapmr.f"
/* L50: */
#line 166 "clapmr.f"
	    ;
#line 166 "clapmr.f"
	}

#line 168 "clapmr.f"
    } else {

/*        Backward permutation */

#line 172 "clapmr.f"
	i__1 = *m;
#line 172 "clapmr.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

#line 174 "clapmr.f"
	    if (k[i__] > 0) {
#line 174 "clapmr.f"
		goto L80;
#line 174 "clapmr.f"
	    }

#line 177 "clapmr.f"
	    k[i__] = -k[i__];
#line 178 "clapmr.f"
	    j = k[i__];
#line 179 "clapmr.f"
L60:
#line 180 "clapmr.f"
	    if (j == i__) {
#line 180 "clapmr.f"
		goto L80;
#line 180 "clapmr.f"
	    }

#line 183 "clapmr.f"
	    i__2 = *n;
#line 183 "clapmr.f"
	    for (jj = 1; jj <= i__2; ++jj) {
#line 184 "clapmr.f"
		i__3 = i__ + jj * x_dim1;
#line 184 "clapmr.f"
		temp.r = x[i__3].r, temp.i = x[i__3].i;
#line 185 "clapmr.f"
		i__3 = i__ + jj * x_dim1;
#line 185 "clapmr.f"
		i__4 = j + jj * x_dim1;
#line 185 "clapmr.f"
		x[i__3].r = x[i__4].r, x[i__3].i = x[i__4].i;
#line 186 "clapmr.f"
		i__3 = j + jj * x_dim1;
#line 186 "clapmr.f"
		x[i__3].r = temp.r, x[i__3].i = temp.i;
#line 187 "clapmr.f"
/* L70: */
#line 187 "clapmr.f"
	    }

#line 189 "clapmr.f"
	    k[j] = -k[j];
#line 190 "clapmr.f"
	    j = k[j];
#line 191 "clapmr.f"
	    goto L60;

#line 193 "clapmr.f"
L80:

#line 195 "clapmr.f"
/* L90: */
#line 195 "clapmr.f"
	    ;
#line 195 "clapmr.f"
	}

#line 197 "clapmr.f"
    }

#line 199 "clapmr.f"
    return 0;

/*     End of ZLAPMT */

} /* clapmr_ */

