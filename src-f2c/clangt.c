#line 1 "clangt.f"
/* clangt.f -- translated by f2c (version 20100827).
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

#line 1 "clangt.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CLANGT returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute 
value of any element of a general tridiagonal matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLANGT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clangt.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clangt.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clangt.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL             FUNCTION CLANGT( NORM, N, DL, D, DU ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          NORM */
/*       INTEGER            N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            D( * ), DL( * ), DU( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLANGT  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the  element of  largest absolute value  of a */
/* > complex tridiagonal matrix A. */
/* > \endverbatim */
/* > */
/* > \return CLANGT */
/* > \verbatim */
/* > */
/* >    CLANGT = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
/* >             ( */
/* >             ( norm1(A),         NORM = '1', 'O' or 'o' */
/* >             ( */
/* >             ( normI(A),         NORM = 'I' or 'i' */
/* >             ( */
/* >             ( normF(A),         NORM = 'F', 'f', 'E' or 'e' */
/* > */
/* > where  norm1  denotes the  one norm of a matrix (maximum column sum), */
/* > normI  denotes the  infinity norm  of a matrix  (maximum row sum) and */
/* > normF  denotes the  Frobenius norm of a matrix (square root of sum of */
/* > squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] NORM */
/* > \verbatim */
/* >          NORM is CHARACTER*1 */
/* >          Specifies the value to be returned in CLANGT as described */
/* >          above. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0.  When N = 0, CLANGT is */
/* >          set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] DL */
/* > \verbatim */
/* >          DL is COMPLEX array, dimension (N-1) */
/* >          The (n-1) sub-diagonal elements of A. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is COMPLEX array, dimension (N) */
/* >          The diagonal elements of A. */
/* > \endverbatim */
/* > */
/* > \param[in] DU */
/* > \verbatim */
/* >          DU is COMPLEX array, dimension (N-1) */
/* >          The (n-1) super-diagonal elements of A. */
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
doublereal clangt_(char *norm, integer *n, doublecomplex *dl, doublecomplex *
	d__, doublecomplex *du, ftnlen norm_len)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double z_abs(doublecomplex *), sqrt(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal sum, temp, scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal anorm;
    extern /* Subroutine */ int classq_(integer *, doublecomplex *, integer *,
	     doublereal *, doublereal *);
    extern logical sisnan_(doublereal *);


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
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 144 "clangt.f"
    /* Parameter adjustments */
#line 144 "clangt.f"
    --du;
#line 144 "clangt.f"
    --d__;
#line 144 "clangt.f"
    --dl;
#line 144 "clangt.f"

#line 144 "clangt.f"
    /* Function Body */
#line 144 "clangt.f"
    if (*n <= 0) {
#line 145 "clangt.f"
	anorm = 0.;
#line 146 "clangt.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

#line 150 "clangt.f"
	anorm = z_abs(&d__[*n]);
#line 151 "clangt.f"
	i__1 = *n - 1;
#line 151 "clangt.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 152 "clangt.f"
	    d__1 = z_abs(&dl[i__]);
#line 152 "clangt.f"
	    if (anorm < z_abs(&dl[i__]) || sisnan_(&d__1)) {
#line 152 "clangt.f"
		anorm = z_abs(&dl[i__]);
#line 152 "clangt.f"
	    }
#line 154 "clangt.f"
	    d__1 = z_abs(&d__[i__]);
#line 154 "clangt.f"
	    if (anorm < z_abs(&d__[i__]) || sisnan_(&d__1)) {
#line 154 "clangt.f"
		anorm = z_abs(&d__[i__]);
#line 154 "clangt.f"
	    }
#line 156 "clangt.f"
	    d__1 = z_abs(&du[i__]);
#line 156 "clangt.f"
	    if (anorm < z_abs(&du[i__]) || sisnan_(&d__1)) {
#line 156 "clangt.f"
		anorm = z_abs(&du[i__]);
#line 156 "clangt.f"
	    }
#line 158 "clangt.f"
/* L10: */
#line 158 "clangt.f"
	}
#line 159 "clangt.f"
    } else if (lsame_(norm, "O", (ftnlen)1, (ftnlen)1) || *(unsigned char *)
	    norm == '1') {

/*        Find norm1(A). */

#line 163 "clangt.f"
	if (*n == 1) {
#line 164 "clangt.f"
	    anorm = z_abs(&d__[1]);
#line 165 "clangt.f"
	} else {
#line 166 "clangt.f"
	    anorm = z_abs(&d__[1]) + z_abs(&dl[1]);
#line 167 "clangt.f"
	    temp = z_abs(&d__[*n]) + z_abs(&du[*n - 1]);
#line 168 "clangt.f"
	    if (anorm < temp || sisnan_(&temp)) {
#line 168 "clangt.f"
		anorm = temp;
#line 168 "clangt.f"
	    }
#line 169 "clangt.f"
	    i__1 = *n - 1;
#line 169 "clangt.f"
	    for (i__ = 2; i__ <= i__1; ++i__) {
#line 170 "clangt.f"
		temp = z_abs(&d__[i__]) + z_abs(&dl[i__]) + z_abs(&du[i__ - 1]
			);
#line 171 "clangt.f"
		if (anorm < temp || sisnan_(&temp)) {
#line 171 "clangt.f"
		    anorm = temp;
#line 171 "clangt.f"
		}
#line 172 "clangt.f"
/* L20: */
#line 172 "clangt.f"
	    }
#line 173 "clangt.f"
	}
#line 174 "clangt.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {

/*        Find normI(A). */

#line 178 "clangt.f"
	if (*n == 1) {
#line 179 "clangt.f"
	    anorm = z_abs(&d__[1]);
#line 180 "clangt.f"
	} else {
#line 181 "clangt.f"
	    anorm = z_abs(&d__[1]) + z_abs(&du[1]);
#line 182 "clangt.f"
	    temp = z_abs(&d__[*n]) + z_abs(&dl[*n - 1]);
#line 183 "clangt.f"
	    if (anorm < temp || sisnan_(&temp)) {
#line 183 "clangt.f"
		anorm = temp;
#line 183 "clangt.f"
	    }
#line 184 "clangt.f"
	    i__1 = *n - 1;
#line 184 "clangt.f"
	    for (i__ = 2; i__ <= i__1; ++i__) {
#line 185 "clangt.f"
		temp = z_abs(&d__[i__]) + z_abs(&du[i__]) + z_abs(&dl[i__ - 1]
			);
#line 186 "clangt.f"
		if (anorm < temp || sisnan_(&temp)) {
#line 186 "clangt.f"
		    anorm = temp;
#line 186 "clangt.f"
		}
#line 187 "clangt.f"
/* L30: */
#line 187 "clangt.f"
	    }
#line 188 "clangt.f"
	}
#line 189 "clangt.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

#line 193 "clangt.f"
	scale = 0.;
#line 194 "clangt.f"
	sum = 1.;
#line 195 "clangt.f"
	classq_(n, &d__[1], &c__1, &scale, &sum);
#line 196 "clangt.f"
	if (*n > 1) {
#line 197 "clangt.f"
	    i__1 = *n - 1;
#line 197 "clangt.f"
	    classq_(&i__1, &dl[1], &c__1, &scale, &sum);
#line 198 "clangt.f"
	    i__1 = *n - 1;
#line 198 "clangt.f"
	    classq_(&i__1, &du[1], &c__1, &scale, &sum);
#line 199 "clangt.f"
	}
#line 200 "clangt.f"
	anorm = scale * sqrt(sum);
#line 201 "clangt.f"
    }

#line 203 "clangt.f"
    ret_val = anorm;
#line 204 "clangt.f"
    return ret_val;

/*     End of CLANGT */

} /* clangt_ */

