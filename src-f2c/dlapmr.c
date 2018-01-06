#line 1 "dlapmr.f"
/* dlapmr.f -- translated by f2c (version 20100827).
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

#line 1 "dlapmr.f"
/* > \brief \b DLAPMR rearranges rows of a matrix as specified by a permutation vector. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLAPMR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlapmr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlapmr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlapmr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLAPMR( FORWRD, M, N, X, LDX, K ) */

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
/* > DLAPMR rearranges the rows of the M by N matrix X as specified */
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

/* > \ingroup doubleOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int dlapmr_(logical *forwrd, integer *m, integer *n, 
	doublereal *x, integer *ldx, integer *k)
{
    /* System generated locals */
    integer x_dim1, x_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, jj, in;
    static doublereal temp;


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

#line 129 "dlapmr.f"
    /* Parameter adjustments */
#line 129 "dlapmr.f"
    x_dim1 = *ldx;
#line 129 "dlapmr.f"
    x_offset = 1 + x_dim1;
#line 129 "dlapmr.f"
    x -= x_offset;
#line 129 "dlapmr.f"
    --k;
#line 129 "dlapmr.f"

#line 129 "dlapmr.f"
    /* Function Body */
#line 129 "dlapmr.f"
    if (*m <= 1) {
#line 129 "dlapmr.f"
	return 0;
#line 129 "dlapmr.f"
    }

#line 132 "dlapmr.f"
    i__1 = *m;
#line 132 "dlapmr.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 133 "dlapmr.f"
	k[i__] = -k[i__];
#line 134 "dlapmr.f"
/* L10: */
#line 134 "dlapmr.f"
    }

#line 136 "dlapmr.f"
    if (*forwrd) {

/*        Forward permutation */

#line 140 "dlapmr.f"
	i__1 = *m;
#line 140 "dlapmr.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

#line 142 "dlapmr.f"
	    if (k[i__] > 0) {
#line 142 "dlapmr.f"
		goto L40;
#line 142 "dlapmr.f"
	    }

#line 145 "dlapmr.f"
	    j = i__;
#line 146 "dlapmr.f"
	    k[j] = -k[j];
#line 147 "dlapmr.f"
	    in = k[j];

#line 149 "dlapmr.f"
L20:
#line 150 "dlapmr.f"
	    if (k[in] > 0) {
#line 150 "dlapmr.f"
		goto L40;
#line 150 "dlapmr.f"
	    }

#line 153 "dlapmr.f"
	    i__2 = *n;
#line 153 "dlapmr.f"
	    for (jj = 1; jj <= i__2; ++jj) {
#line 154 "dlapmr.f"
		temp = x[j + jj * x_dim1];
#line 155 "dlapmr.f"
		x[j + jj * x_dim1] = x[in + jj * x_dim1];
#line 156 "dlapmr.f"
		x[in + jj * x_dim1] = temp;
#line 157 "dlapmr.f"
/* L30: */
#line 157 "dlapmr.f"
	    }

#line 159 "dlapmr.f"
	    k[in] = -k[in];
#line 160 "dlapmr.f"
	    j = in;
#line 161 "dlapmr.f"
	    in = k[in];
#line 162 "dlapmr.f"
	    goto L20;

#line 164 "dlapmr.f"
L40:

#line 166 "dlapmr.f"
/* L50: */
#line 166 "dlapmr.f"
	    ;
#line 166 "dlapmr.f"
	}

#line 168 "dlapmr.f"
    } else {

/*        Backward permutation */

#line 172 "dlapmr.f"
	i__1 = *m;
#line 172 "dlapmr.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

#line 174 "dlapmr.f"
	    if (k[i__] > 0) {
#line 174 "dlapmr.f"
		goto L80;
#line 174 "dlapmr.f"
	    }

#line 177 "dlapmr.f"
	    k[i__] = -k[i__];
#line 178 "dlapmr.f"
	    j = k[i__];
#line 179 "dlapmr.f"
L60:
#line 180 "dlapmr.f"
	    if (j == i__) {
#line 180 "dlapmr.f"
		goto L80;
#line 180 "dlapmr.f"
	    }

#line 183 "dlapmr.f"
	    i__2 = *n;
#line 183 "dlapmr.f"
	    for (jj = 1; jj <= i__2; ++jj) {
#line 184 "dlapmr.f"
		temp = x[i__ + jj * x_dim1];
#line 185 "dlapmr.f"
		x[i__ + jj * x_dim1] = x[j + jj * x_dim1];
#line 186 "dlapmr.f"
		x[j + jj * x_dim1] = temp;
#line 187 "dlapmr.f"
/* L70: */
#line 187 "dlapmr.f"
	    }

#line 189 "dlapmr.f"
	    k[j] = -k[j];
#line 190 "dlapmr.f"
	    j = k[j];
#line 191 "dlapmr.f"
	    goto L60;

#line 193 "dlapmr.f"
L80:

#line 195 "dlapmr.f"
/* L90: */
#line 195 "dlapmr.f"
	    ;
#line 195 "dlapmr.f"
	}

#line 197 "dlapmr.f"
    }

#line 199 "dlapmr.f"
    return 0;

/*     End of ZLAPMT */

} /* dlapmr_ */

