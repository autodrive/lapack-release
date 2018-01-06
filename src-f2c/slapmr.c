#line 1 "slapmr.f"
/* slapmr.f -- translated by f2c (version 20100827).
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

#line 1 "slapmr.f"
/* > \brief \b SLAPMR rearranges rows of a matrix as specified by a permutation vector. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLAPMR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slapmr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slapmr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slapmr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLAPMR( FORWRD, M, N, X, LDX, K ) */

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
/* > SLAPMR rearranges the rows of the M by N matrix X as specified */
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

/* > \ingroup realOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int slapmr_(logical *forwrd, integer *m, integer *n, 
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

#line 129 "slapmr.f"
    /* Parameter adjustments */
#line 129 "slapmr.f"
    x_dim1 = *ldx;
#line 129 "slapmr.f"
    x_offset = 1 + x_dim1;
#line 129 "slapmr.f"
    x -= x_offset;
#line 129 "slapmr.f"
    --k;
#line 129 "slapmr.f"

#line 129 "slapmr.f"
    /* Function Body */
#line 129 "slapmr.f"
    if (*m <= 1) {
#line 129 "slapmr.f"
	return 0;
#line 129 "slapmr.f"
    }

#line 132 "slapmr.f"
    i__1 = *m;
#line 132 "slapmr.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 133 "slapmr.f"
	k[i__] = -k[i__];
#line 134 "slapmr.f"
/* L10: */
#line 134 "slapmr.f"
    }

#line 136 "slapmr.f"
    if (*forwrd) {

/*        Forward permutation */

#line 140 "slapmr.f"
	i__1 = *m;
#line 140 "slapmr.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

#line 142 "slapmr.f"
	    if (k[i__] > 0) {
#line 142 "slapmr.f"
		goto L40;
#line 142 "slapmr.f"
	    }

#line 145 "slapmr.f"
	    j = i__;
#line 146 "slapmr.f"
	    k[j] = -k[j];
#line 147 "slapmr.f"
	    in = k[j];

#line 149 "slapmr.f"
L20:
#line 150 "slapmr.f"
	    if (k[in] > 0) {
#line 150 "slapmr.f"
		goto L40;
#line 150 "slapmr.f"
	    }

#line 153 "slapmr.f"
	    i__2 = *n;
#line 153 "slapmr.f"
	    for (jj = 1; jj <= i__2; ++jj) {
#line 154 "slapmr.f"
		temp = x[j + jj * x_dim1];
#line 155 "slapmr.f"
		x[j + jj * x_dim1] = x[in + jj * x_dim1];
#line 156 "slapmr.f"
		x[in + jj * x_dim1] = temp;
#line 157 "slapmr.f"
/* L30: */
#line 157 "slapmr.f"
	    }

#line 159 "slapmr.f"
	    k[in] = -k[in];
#line 160 "slapmr.f"
	    j = in;
#line 161 "slapmr.f"
	    in = k[in];
#line 162 "slapmr.f"
	    goto L20;

#line 164 "slapmr.f"
L40:

#line 166 "slapmr.f"
/* L50: */
#line 166 "slapmr.f"
	    ;
#line 166 "slapmr.f"
	}

#line 168 "slapmr.f"
    } else {

/*        Backward permutation */

#line 172 "slapmr.f"
	i__1 = *m;
#line 172 "slapmr.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

#line 174 "slapmr.f"
	    if (k[i__] > 0) {
#line 174 "slapmr.f"
		goto L80;
#line 174 "slapmr.f"
	    }

#line 177 "slapmr.f"
	    k[i__] = -k[i__];
#line 178 "slapmr.f"
	    j = k[i__];
#line 179 "slapmr.f"
L60:
#line 180 "slapmr.f"
	    if (j == i__) {
#line 180 "slapmr.f"
		goto L80;
#line 180 "slapmr.f"
	    }

#line 183 "slapmr.f"
	    i__2 = *n;
#line 183 "slapmr.f"
	    for (jj = 1; jj <= i__2; ++jj) {
#line 184 "slapmr.f"
		temp = x[i__ + jj * x_dim1];
#line 185 "slapmr.f"
		x[i__ + jj * x_dim1] = x[j + jj * x_dim1];
#line 186 "slapmr.f"
		x[j + jj * x_dim1] = temp;
#line 187 "slapmr.f"
/* L70: */
#line 187 "slapmr.f"
	    }

#line 189 "slapmr.f"
	    k[j] = -k[j];
#line 190 "slapmr.f"
	    j = k[j];
#line 191 "slapmr.f"
	    goto L60;

#line 193 "slapmr.f"
L80:

#line 195 "slapmr.f"
/* L90: */
#line 195 "slapmr.f"
	    ;
#line 195 "slapmr.f"
	}

#line 197 "slapmr.f"
    }

#line 199 "slapmr.f"
    return 0;

/*     End of ZLAPMT */

} /* slapmr_ */

