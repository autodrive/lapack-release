#line 1 "cpptri.f"
/* cpptri.f -- translated by f2c (version 20100827).
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

#line 1 "cpptri.f"
/* Table of constant values */

static doublereal c_b8 = 1.;
static integer c__1 = 1;

/* > \brief \b CPPTRI */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CPPTRI + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpptri.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpptri.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpptri.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CPPTRI( UPLO, N, AP, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            AP( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CPPTRI computes the inverse of a complex Hermitian positive definite */
/* > matrix A using the Cholesky factorization A = U**H*U or A = L*L**H */
/* > computed by CPPTRF. */
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
/* >          AP is COMPLEX array, dimension (N*(N+1)/2) */
/* >          On entry, the triangular factor U or L from the Cholesky */
/* >          factorization A = U**H*U or A = L*L**H, packed columnwise as */
/* >          a linear array.  The j-th column of U or L is stored in the */
/* >          array AP as follows: */
/* >          if UPLO = 'U', AP(i + (j-1)*j/2) = U(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = L(i,j) for j<=i<=n. */
/* > */
/* >          On exit, the upper or lower triangle of the (Hermitian) */
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

/* > \ingroup complexOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int cpptri_(char *uplo, integer *n, doublecomplex *ap, 
	integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;
    doublecomplex z__1;

    /* Local variables */
    static integer j, jc, jj;
    static doublereal ajj;
    static integer jjn;
    extern /* Subroutine */ int chpr_(char *, integer *, doublereal *, 
	    doublecomplex *, integer *, doublecomplex *, ftnlen);
    extern /* Double Complex */ VOID cdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int ctpmv_(char *, char *, char *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen, ftnlen, 
	    ftnlen);
    static logical upper;
    extern /* Subroutine */ int csscal_(integer *, doublereal *, 
	    doublecomplex *, integer *), xerbla_(char *, integer *, ftnlen), 
	    ctptri_(char *, char *, integer *, doublecomplex *, integer *, 
	    ftnlen, ftnlen);


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
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 135 "cpptri.f"
    /* Parameter adjustments */
#line 135 "cpptri.f"
    --ap;
#line 135 "cpptri.f"

#line 135 "cpptri.f"
    /* Function Body */
#line 135 "cpptri.f"
    *info = 0;
#line 136 "cpptri.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 137 "cpptri.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 138 "cpptri.f"
	*info = -1;
#line 139 "cpptri.f"
    } else if (*n < 0) {
#line 140 "cpptri.f"
	*info = -2;
#line 141 "cpptri.f"
    }
#line 142 "cpptri.f"
    if (*info != 0) {
#line 143 "cpptri.f"
	i__1 = -(*info);
#line 143 "cpptri.f"
	xerbla_("CPPTRI", &i__1, (ftnlen)6);
#line 144 "cpptri.f"
	return 0;
#line 145 "cpptri.f"
    }

/*     Quick return if possible */

#line 149 "cpptri.f"
    if (*n == 0) {
#line 149 "cpptri.f"
	return 0;
#line 149 "cpptri.f"
    }

/*     Invert the triangular Cholesky factor U or L. */

#line 154 "cpptri.f"
    ctptri_(uplo, "Non-unit", n, &ap[1], info, (ftnlen)1, (ftnlen)8);
#line 155 "cpptri.f"
    if (*info > 0) {
#line 155 "cpptri.f"
	return 0;
#line 155 "cpptri.f"
    }
#line 157 "cpptri.f"
    if (upper) {

/*        Compute the product inv(U) * inv(U)**H. */

#line 161 "cpptri.f"
	jj = 0;
#line 162 "cpptri.f"
	i__1 = *n;
#line 162 "cpptri.f"
	for (j = 1; j <= i__1; ++j) {
#line 163 "cpptri.f"
	    jc = jj + 1;
#line 164 "cpptri.f"
	    jj += j;
#line 165 "cpptri.f"
	    if (j > 1) {
#line 165 "cpptri.f"
		i__2 = j - 1;
#line 165 "cpptri.f"
		chpr_("Upper", &i__2, &c_b8, &ap[jc], &c__1, &ap[1], (ftnlen)
			5);
#line 165 "cpptri.f"
	    }
#line 167 "cpptri.f"
	    i__2 = jj;
#line 167 "cpptri.f"
	    ajj = ap[i__2].r;
#line 168 "cpptri.f"
	    csscal_(&j, &ajj, &ap[jc], &c__1);
#line 169 "cpptri.f"
/* L10: */
#line 169 "cpptri.f"
	}

#line 171 "cpptri.f"
    } else {

/*        Compute the product inv(L)**H * inv(L). */

#line 175 "cpptri.f"
	jj = 1;
#line 176 "cpptri.f"
	i__1 = *n;
#line 176 "cpptri.f"
	for (j = 1; j <= i__1; ++j) {
#line 177 "cpptri.f"
	    jjn = jj + *n - j + 1;
#line 178 "cpptri.f"
	    i__2 = jj;
#line 178 "cpptri.f"
	    i__3 = *n - j + 1;
#line 178 "cpptri.f"
	    cdotc_(&z__1, &i__3, &ap[jj], &c__1, &ap[jj], &c__1);
#line 178 "cpptri.f"
	    d__1 = z__1.r;
#line 178 "cpptri.f"
	    ap[i__2].r = d__1, ap[i__2].i = 0.;
#line 179 "cpptri.f"
	    if (j < *n) {
#line 179 "cpptri.f"
		i__2 = *n - j;
#line 179 "cpptri.f"
		ctpmv_("Lower", "Conjugate transpose", "Non-unit", &i__2, &ap[
			jjn], &ap[jj + 1], &c__1, (ftnlen)5, (ftnlen)19, (
			ftnlen)8);
#line 179 "cpptri.f"
	    }
#line 182 "cpptri.f"
	    jj = jjn;
#line 183 "cpptri.f"
/* L20: */
#line 183 "cpptri.f"
	}
#line 184 "cpptri.f"
    }

#line 186 "cpptri.f"
    return 0;

/*     End of CPPTRI */

} /* cpptri_ */

