#line 1 "sgetri.f"
/* sgetri.f -- translated by f2c (version 20100827).
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

#line 1 "sgetri.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static doublereal c_b20 = -1.;
static doublereal c_b22 = 1.;

/* > \brief \b SGETRI */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGETRI + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgetri.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgetri.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgetri.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       REAL               A( LDA, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGETRI computes the inverse of a matrix using the LU factorization */
/* > computed by SGETRF. */
/* > */
/* > This method inverts U and then computes inv(A) by solving the system */
/* > inv(A)*L = inv(U) for inv(A). */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA,N) */
/* >          On entry, the factors L and U from the factorization */
/* >          A = P*L*U as computed by SGETRF. */
/* >          On exit, if INFO = 0, the inverse of the original matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          The pivot indices from SGETRF; for 1<=i<=N, row i of the */
/* >          matrix was interchanged with row IPIV(i). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO=0, then WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK.  LWORK >= max(1,N). */
/* >          For optimal performance LWORK >= N*NB, where NB is */
/* >          the optimal blocksize returned by ILAENV. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is */
/* >                singular and its inverse could not be computed. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup realGEcomputational */

/*  ===================================================================== */
/* Subroutine */ int sgetri_(integer *n, doublereal *a, integer *lda, integer 
	*ipiv, doublereal *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, jb, nb, jj, jp, nn, iws, nbmin;
    extern /* Subroutine */ int sgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen),
	     sgemv_(char *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, ftnlen), sswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), strsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), xerbla_(
	    char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer ldwork, lwkopt;
    static logical lquery;
    extern /* Subroutine */ int strtri_(char *, char *, integer *, doublereal 
	    *, integer *, integer *, ftnlen, ftnlen);


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

#line 155 "sgetri.f"
    /* Parameter adjustments */
#line 155 "sgetri.f"
    a_dim1 = *lda;
#line 155 "sgetri.f"
    a_offset = 1 + a_dim1;
#line 155 "sgetri.f"
    a -= a_offset;
#line 155 "sgetri.f"
    --ipiv;
#line 155 "sgetri.f"
    --work;
#line 155 "sgetri.f"

#line 155 "sgetri.f"
    /* Function Body */
#line 155 "sgetri.f"
    *info = 0;
#line 156 "sgetri.f"
    nb = ilaenv_(&c__1, "SGETRI", " ", n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)1);
#line 157 "sgetri.f"
    lwkopt = *n * nb;
#line 158 "sgetri.f"
    work[1] = (doublereal) lwkopt;
#line 159 "sgetri.f"
    lquery = *lwork == -1;
#line 160 "sgetri.f"
    if (*n < 0) {
#line 161 "sgetri.f"
	*info = -1;
#line 162 "sgetri.f"
    } else if (*lda < max(1,*n)) {
#line 163 "sgetri.f"
	*info = -3;
#line 164 "sgetri.f"
    } else if (*lwork < max(1,*n) && ! lquery) {
#line 165 "sgetri.f"
	*info = -6;
#line 166 "sgetri.f"
    }
#line 167 "sgetri.f"
    if (*info != 0) {
#line 168 "sgetri.f"
	i__1 = -(*info);
#line 168 "sgetri.f"
	xerbla_("SGETRI", &i__1, (ftnlen)6);
#line 169 "sgetri.f"
	return 0;
#line 170 "sgetri.f"
    } else if (lquery) {
#line 171 "sgetri.f"
	return 0;
#line 172 "sgetri.f"
    }

/*     Quick return if possible */

#line 176 "sgetri.f"
    if (*n == 0) {
#line 176 "sgetri.f"
	return 0;
#line 176 "sgetri.f"
    }

/*     Form inv(U).  If INFO > 0 from STRTRI, then U is singular, */
/*     and the inverse is not computed. */

#line 182 "sgetri.f"
    strtri_("Upper", "Non-unit", n, &a[a_offset], lda, info, (ftnlen)5, (
	    ftnlen)8);
#line 183 "sgetri.f"
    if (*info > 0) {
#line 183 "sgetri.f"
	return 0;
#line 183 "sgetri.f"
    }

#line 186 "sgetri.f"
    nbmin = 2;
#line 187 "sgetri.f"
    ldwork = *n;
#line 188 "sgetri.f"
    if (nb > 1 && nb < *n) {
/* Computing MAX */
#line 189 "sgetri.f"
	i__1 = ldwork * nb;
#line 189 "sgetri.f"
	iws = max(i__1,1);
#line 190 "sgetri.f"
	if (*lwork < iws) {
#line 191 "sgetri.f"
	    nb = *lwork / ldwork;
/* Computing MAX */
#line 192 "sgetri.f"
	    i__1 = 2, i__2 = ilaenv_(&c__2, "SGETRI", " ", n, &c_n1, &c_n1, &
		    c_n1, (ftnlen)6, (ftnlen)1);
#line 192 "sgetri.f"
	    nbmin = max(i__1,i__2);
#line 193 "sgetri.f"
	}
#line 194 "sgetri.f"
    } else {
#line 195 "sgetri.f"
	iws = *n;
#line 196 "sgetri.f"
    }

/*     Solve the equation inv(A)*L = inv(U) for inv(A). */

#line 200 "sgetri.f"
    if (nb < nbmin || nb >= *n) {

/*        Use unblocked code. */

#line 204 "sgetri.f"
	for (j = *n; j >= 1; --j) {

/*           Copy current column of L to WORK and replace with zeros. */

#line 208 "sgetri.f"
	    i__1 = *n;
#line 208 "sgetri.f"
	    for (i__ = j + 1; i__ <= i__1; ++i__) {
#line 209 "sgetri.f"
		work[i__] = a[i__ + j * a_dim1];
#line 210 "sgetri.f"
		a[i__ + j * a_dim1] = 0.;
#line 211 "sgetri.f"
/* L10: */
#line 211 "sgetri.f"
	    }

/*           Compute current column of inv(A). */

#line 215 "sgetri.f"
	    if (j < *n) {
#line 215 "sgetri.f"
		i__1 = *n - j;
#line 215 "sgetri.f"
		sgemv_("No transpose", n, &i__1, &c_b20, &a[(j + 1) * a_dim1 
			+ 1], lda, &work[j + 1], &c__1, &c_b22, &a[j * a_dim1 
			+ 1], &c__1, (ftnlen)12);
#line 215 "sgetri.f"
	    }
#line 218 "sgetri.f"
/* L20: */
#line 218 "sgetri.f"
	}
#line 219 "sgetri.f"
    } else {

/*        Use blocked code. */

#line 223 "sgetri.f"
	nn = (*n - 1) / nb * nb + 1;
#line 224 "sgetri.f"
	i__1 = -nb;
#line 224 "sgetri.f"
	for (j = nn; i__1 < 0 ? j >= 1 : j <= 1; j += i__1) {
/* Computing MIN */
#line 225 "sgetri.f"
	    i__2 = nb, i__3 = *n - j + 1;
#line 225 "sgetri.f"
	    jb = min(i__2,i__3);

/*           Copy current block column of L to WORK and replace with */
/*           zeros. */

#line 230 "sgetri.f"
	    i__2 = j + jb - 1;
#line 230 "sgetri.f"
	    for (jj = j; jj <= i__2; ++jj) {
#line 231 "sgetri.f"
		i__3 = *n;
#line 231 "sgetri.f"
		for (i__ = jj + 1; i__ <= i__3; ++i__) {
#line 232 "sgetri.f"
		    work[i__ + (jj - j) * ldwork] = a[i__ + jj * a_dim1];
#line 233 "sgetri.f"
		    a[i__ + jj * a_dim1] = 0.;
#line 234 "sgetri.f"
/* L30: */
#line 234 "sgetri.f"
		}
#line 235 "sgetri.f"
/* L40: */
#line 235 "sgetri.f"
	    }

/*           Compute current block column of inv(A). */

#line 239 "sgetri.f"
	    if (j + jb <= *n) {
#line 239 "sgetri.f"
		i__2 = *n - j - jb + 1;
#line 239 "sgetri.f"
		sgemm_("No transpose", "No transpose", n, &jb, &i__2, &c_b20, 
			&a[(j + jb) * a_dim1 + 1], lda, &work[j + jb], &
			ldwork, &c_b22, &a[j * a_dim1 + 1], lda, (ftnlen)12, (
			ftnlen)12);
#line 239 "sgetri.f"
	    }
#line 243 "sgetri.f"
	    strsm_("Right", "Lower", "No transpose", "Unit", n, &jb, &c_b22, &
		    work[j], &ldwork, &a[j * a_dim1 + 1], lda, (ftnlen)5, (
		    ftnlen)5, (ftnlen)12, (ftnlen)4);
#line 245 "sgetri.f"
/* L50: */
#line 245 "sgetri.f"
	}
#line 246 "sgetri.f"
    }

/*     Apply column interchanges. */

#line 250 "sgetri.f"
    for (j = *n - 1; j >= 1; --j) {
#line 251 "sgetri.f"
	jp = ipiv[j];
#line 252 "sgetri.f"
	if (jp != j) {
#line 252 "sgetri.f"
	    sswap_(n, &a[j * a_dim1 + 1], &c__1, &a[jp * a_dim1 + 1], &c__1);
#line 252 "sgetri.f"
	}
#line 254 "sgetri.f"
/* L60: */
#line 254 "sgetri.f"
    }

#line 256 "sgetri.f"
    work[1] = (doublereal) iws;
#line 257 "sgetri.f"
    return 0;

/*     End of SGETRI */

} /* sgetri_ */

