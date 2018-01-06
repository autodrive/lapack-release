#line 1 "spptri.f"
/* spptri.f -- translated by f2c (version 20100827).
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

#line 1 "spptri.f"
/* Table of constant values */

static doublereal c_b8 = 1.;
static integer c__1 = 1;

/* > \brief \b SPPTRI */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SPPTRI + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/spptri.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/spptri.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/spptri.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SPPTRI( UPLO, N, AP, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               AP( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SPPTRI computes the inverse of a real symmetric positive definite */
/* > matrix A using the Cholesky factorization A = U**T*U or A = L*L**T */
/* > computed by SPPTRF. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U':  Upper triangular factor is stored in AP; */
/* >          = 'L':  Lower triangular factor is stored in AP. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AP */
/* > \verbatim */
/* >          AP is REAL array, dimension (N*(N+1)/2) */
/* >          On entry, the triangular factor U or L from the Cholesky */
/* >          factorization A = U**T*U or A = L*L**T, packed columnwise as */
/* >          a linear array.  The j-th column of U or L is stored in the */
/* >          array AP as follows: */
/* >          if UPLO = 'U', AP(i + (j-1)*j/2) = U(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = L(i,j) for j<=i<=n. */
/* > */
/* >          On exit, the upper or lower triangle of the (symmetric) */
/* >          inverse of A, overwriting the input factor U or L. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  if INFO = i, the (i,i) element of the factor U or L is */
/* >                zero, and the inverse could not be computed. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup realOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int spptri_(char *uplo, integer *n, doublereal *ap, integer *
	info, ftnlen uplo_len)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer j, jc, jj;
    static doublereal ajj;
    static integer jjn;
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern /* Subroutine */ int sspr_(char *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static logical upper;
    extern /* Subroutine */ int stpmv_(char *, char *, char *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen, ftnlen), 
	    xerbla_(char *, integer *, ftnlen), stptri_(char *, char *, 
	    integer *, doublereal *, integer *, ftnlen, ftnlen);


/*  -- LAPACK computational routine (version 3.7.0) -- */
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
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 132 "spptri.f"
    /* Parameter adjustments */
#line 132 "spptri.f"
    --ap;
#line 132 "spptri.f"

#line 132 "spptri.f"
    /* Function Body */
#line 132 "spptri.f"
    *info = 0;
#line 133 "spptri.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 134 "spptri.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 135 "spptri.f"
	*info = -1;
#line 136 "spptri.f"
    } else if (*n < 0) {
#line 137 "spptri.f"
	*info = -2;
#line 138 "spptri.f"
    }
#line 139 "spptri.f"
    if (*info != 0) {
#line 140 "spptri.f"
	i__1 = -(*info);
#line 140 "spptri.f"
	xerbla_("SPPTRI", &i__1, (ftnlen)6);
#line 141 "spptri.f"
	return 0;
#line 142 "spptri.f"
    }

/*     Quick return if possible */

#line 146 "spptri.f"
    if (*n == 0) {
#line 146 "spptri.f"
	return 0;
#line 146 "spptri.f"
    }

/*     Invert the triangular Cholesky factor U or L. */

#line 151 "spptri.f"
    stptri_(uplo, "Non-unit", n, &ap[1], info, (ftnlen)1, (ftnlen)8);
#line 152 "spptri.f"
    if (*info > 0) {
#line 152 "spptri.f"
	return 0;
#line 152 "spptri.f"
    }

#line 155 "spptri.f"
    if (upper) {

/*        Compute the product inv(U) * inv(U)**T. */

#line 159 "spptri.f"
	jj = 0;
#line 160 "spptri.f"
	i__1 = *n;
#line 160 "spptri.f"
	for (j = 1; j <= i__1; ++j) {
#line 161 "spptri.f"
	    jc = jj + 1;
#line 162 "spptri.f"
	    jj += j;
#line 163 "spptri.f"
	    if (j > 1) {
#line 163 "spptri.f"
		i__2 = j - 1;
#line 163 "spptri.f"
		sspr_("Upper", &i__2, &c_b8, &ap[jc], &c__1, &ap[1], (ftnlen)
			5);
#line 163 "spptri.f"
	    }
#line 165 "spptri.f"
	    ajj = ap[jj];
#line 166 "spptri.f"
	    sscal_(&j, &ajj, &ap[jc], &c__1);
#line 167 "spptri.f"
/* L10: */
#line 167 "spptri.f"
	}

#line 169 "spptri.f"
    } else {

/*        Compute the product inv(L)**T * inv(L). */

#line 173 "spptri.f"
	jj = 1;
#line 174 "spptri.f"
	i__1 = *n;
#line 174 "spptri.f"
	for (j = 1; j <= i__1; ++j) {
#line 175 "spptri.f"
	    jjn = jj + *n - j + 1;
#line 176 "spptri.f"
	    i__2 = *n - j + 1;
#line 176 "spptri.f"
	    ap[jj] = sdot_(&i__2, &ap[jj], &c__1, &ap[jj], &c__1);
#line 177 "spptri.f"
	    if (j < *n) {
#line 177 "spptri.f"
		i__2 = *n - j;
#line 177 "spptri.f"
		stpmv_("Lower", "Transpose", "Non-unit", &i__2, &ap[jjn], &ap[
			jj + 1], &c__1, (ftnlen)5, (ftnlen)9, (ftnlen)8);
#line 177 "spptri.f"
	    }
#line 180 "spptri.f"
	    jj = jjn;
#line 181 "spptri.f"
/* L20: */
#line 181 "spptri.f"
	}
#line 182 "spptri.f"
    }

#line 184 "spptri.f"
    return 0;

/*     End of SPPTRI */

} /* spptri_ */

