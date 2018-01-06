#line 1 "zpptri.f"
/* zpptri.f -- translated by f2c (version 20100827).
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

#line 1 "zpptri.f"
/* Table of constant values */

static doublereal c_b8 = 1.;
static integer c__1 = 1;

/* > \brief \b ZPPTRI */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZPPTRI + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zpptri.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zpptri.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zpptri.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZPPTRI( UPLO, N, AP, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         AP( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZPPTRI computes the inverse of a complex Hermitian positive definite */
/* > matrix A using the Cholesky factorization A = U**H*U or A = L*L**H */
/* > computed by ZPPTRF. */
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
/* >          AP is COMPLEX*16 array, dimension (N*(N+1)/2) */
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

/* > \date November 2011 */

/* > \ingroup complex16OTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int zpptri_(char *uplo, integer *n, doublecomplex *ap, 
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
    extern /* Subroutine */ int zhpr_(char *, integer *, doublereal *, 
	    doublecomplex *, integer *, doublecomplex *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Double Complex */ VOID zdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static logical upper;
    extern /* Subroutine */ int ztpmv_(char *, char *, char *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen, ftnlen, 
	    ftnlen), xerbla_(char *, integer *, ftnlen), zdscal_(integer *, 
	    doublereal *, doublecomplex *, integer *), ztptri_(char *, char *,
	     integer *, doublecomplex *, integer *, ftnlen, ftnlen);


/*  -- LAPACK computational routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

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

#line 135 "zpptri.f"
    /* Parameter adjustments */
#line 135 "zpptri.f"
    --ap;
#line 135 "zpptri.f"

#line 135 "zpptri.f"
    /* Function Body */
#line 135 "zpptri.f"
    *info = 0;
#line 136 "zpptri.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 137 "zpptri.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 138 "zpptri.f"
	*info = -1;
#line 139 "zpptri.f"
    } else if (*n < 0) {
#line 140 "zpptri.f"
	*info = -2;
#line 141 "zpptri.f"
    }
#line 142 "zpptri.f"
    if (*info != 0) {
#line 143 "zpptri.f"
	i__1 = -(*info);
#line 143 "zpptri.f"
	xerbla_("ZPPTRI", &i__1, (ftnlen)6);
#line 144 "zpptri.f"
	return 0;
#line 145 "zpptri.f"
    }

/*     Quick return if possible */

#line 149 "zpptri.f"
    if (*n == 0) {
#line 149 "zpptri.f"
	return 0;
#line 149 "zpptri.f"
    }

/*     Invert the triangular Cholesky factor U or L. */

#line 154 "zpptri.f"
    ztptri_(uplo, "Non-unit", n, &ap[1], info, (ftnlen)1, (ftnlen)8);
#line 155 "zpptri.f"
    if (*info > 0) {
#line 155 "zpptri.f"
	return 0;
#line 155 "zpptri.f"
    }
#line 157 "zpptri.f"
    if (upper) {

/*        Compute the product inv(U) * inv(U)**H. */

#line 161 "zpptri.f"
	jj = 0;
#line 162 "zpptri.f"
	i__1 = *n;
#line 162 "zpptri.f"
	for (j = 1; j <= i__1; ++j) {
#line 163 "zpptri.f"
	    jc = jj + 1;
#line 164 "zpptri.f"
	    jj += j;
#line 165 "zpptri.f"
	    if (j > 1) {
#line 165 "zpptri.f"
		i__2 = j - 1;
#line 165 "zpptri.f"
		zhpr_("Upper", &i__2, &c_b8, &ap[jc], &c__1, &ap[1], (ftnlen)
			5);
#line 165 "zpptri.f"
	    }
#line 167 "zpptri.f"
	    i__2 = jj;
#line 167 "zpptri.f"
	    ajj = ap[i__2].r;
#line 168 "zpptri.f"
	    zdscal_(&j, &ajj, &ap[jc], &c__1);
#line 169 "zpptri.f"
/* L10: */
#line 169 "zpptri.f"
	}

#line 171 "zpptri.f"
    } else {

/*        Compute the product inv(L)**H * inv(L). */

#line 175 "zpptri.f"
	jj = 1;
#line 176 "zpptri.f"
	i__1 = *n;
#line 176 "zpptri.f"
	for (j = 1; j <= i__1; ++j) {
#line 177 "zpptri.f"
	    jjn = jj + *n - j + 1;
#line 178 "zpptri.f"
	    i__2 = jj;
#line 178 "zpptri.f"
	    i__3 = *n - j + 1;
#line 178 "zpptri.f"
	    zdotc_(&z__1, &i__3, &ap[jj], &c__1, &ap[jj], &c__1);
#line 178 "zpptri.f"
	    d__1 = z__1.r;
#line 178 "zpptri.f"
	    ap[i__2].r = d__1, ap[i__2].i = 0.;
#line 179 "zpptri.f"
	    if (j < *n) {
#line 179 "zpptri.f"
		i__2 = *n - j;
#line 179 "zpptri.f"
		ztpmv_("Lower", "Conjugate transpose", "Non-unit", &i__2, &ap[
			jjn], &ap[jj + 1], &c__1, (ftnlen)5, (ftnlen)19, (
			ftnlen)8);
#line 179 "zpptri.f"
	    }
#line 182 "zpptri.f"
	    jj = jjn;
#line 183 "zpptri.f"
/* L20: */
#line 183 "zpptri.f"
	}
#line 184 "zpptri.f"
    }

#line 186 "zpptri.f"
    return 0;

/*     End of ZPPTRI */

} /* zpptri_ */

