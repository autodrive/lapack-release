#line 1 "cgetri.f"
/* cgetri.f -- translated by f2c (version 20100827).
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

#line 1 "cgetri.f"
/* Table of constant values */

static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;

/* > \brief \b CGETRI */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGETRI + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgetri.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgetri.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgetri.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX            A( LDA, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGETRI computes the inverse of a matrix using the LU factorization */
/* > computed by CGETRF. */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          On entry, the factors L and U from the factorization */
/* >          A = P*L*U as computed by CGETRF. */
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
/* >          The pivot indices from CGETRF; for 1<=i<=N, row i of the */
/* >          matrix was interchanged with row IPIV(i). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
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

/* > \ingroup complexGEcomputational */

/*  ===================================================================== */
/* Subroutine */ int cgetri_(integer *n, doublecomplex *a, integer *lda, 
	integer *ipiv, doublecomplex *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1;

    /* Local variables */
    static integer i__, j, jb, nb, jj, jp, nn, iws;
    extern /* Subroutine */ int cgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen), cgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen);
    static integer nbmin;
    extern /* Subroutine */ int cswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), ctrsm_(char *, char *, char *, char *
	    , integer *, integer *, doublecomplex *, doublecomplex *, integer 
	    *, doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer ldwork;
    extern /* Subroutine */ int ctrtri_(char *, char *, integer *, 
	    doublecomplex *, integer *, integer *, ftnlen, ftnlen);
    static integer lwkopt;
    static logical lquery;


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

#line 156 "cgetri.f"
    /* Parameter adjustments */
#line 156 "cgetri.f"
    a_dim1 = *lda;
#line 156 "cgetri.f"
    a_offset = 1 + a_dim1;
#line 156 "cgetri.f"
    a -= a_offset;
#line 156 "cgetri.f"
    --ipiv;
#line 156 "cgetri.f"
    --work;
#line 156 "cgetri.f"

#line 156 "cgetri.f"
    /* Function Body */
#line 156 "cgetri.f"
    *info = 0;
#line 157 "cgetri.f"
    nb = ilaenv_(&c__1, "CGETRI", " ", n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)1);
#line 158 "cgetri.f"
    lwkopt = *n * nb;
#line 159 "cgetri.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 160 "cgetri.f"
    lquery = *lwork == -1;
#line 161 "cgetri.f"
    if (*n < 0) {
#line 162 "cgetri.f"
	*info = -1;
#line 163 "cgetri.f"
    } else if (*lda < max(1,*n)) {
#line 164 "cgetri.f"
	*info = -3;
#line 165 "cgetri.f"
    } else if (*lwork < max(1,*n) && ! lquery) {
#line 166 "cgetri.f"
	*info = -6;
#line 167 "cgetri.f"
    }
#line 168 "cgetri.f"
    if (*info != 0) {
#line 169 "cgetri.f"
	i__1 = -(*info);
#line 169 "cgetri.f"
	xerbla_("CGETRI", &i__1, (ftnlen)6);
#line 170 "cgetri.f"
	return 0;
#line 171 "cgetri.f"
    } else if (lquery) {
#line 172 "cgetri.f"
	return 0;
#line 173 "cgetri.f"
    }

/*     Quick return if possible */

#line 177 "cgetri.f"
    if (*n == 0) {
#line 177 "cgetri.f"
	return 0;
#line 177 "cgetri.f"
    }

/*     Form inv(U).  If INFO > 0 from CTRTRI, then U is singular, */
/*     and the inverse is not computed. */

#line 183 "cgetri.f"
    ctrtri_("Upper", "Non-unit", n, &a[a_offset], lda, info, (ftnlen)5, (
	    ftnlen)8);
#line 184 "cgetri.f"
    if (*info > 0) {
#line 184 "cgetri.f"
	return 0;
#line 184 "cgetri.f"
    }

#line 187 "cgetri.f"
    nbmin = 2;
#line 188 "cgetri.f"
    ldwork = *n;
#line 189 "cgetri.f"
    if (nb > 1 && nb < *n) {
/* Computing MAX */
#line 190 "cgetri.f"
	i__1 = ldwork * nb;
#line 190 "cgetri.f"
	iws = max(i__1,1);
#line 191 "cgetri.f"
	if (*lwork < iws) {
#line 192 "cgetri.f"
	    nb = *lwork / ldwork;
/* Computing MAX */
#line 193 "cgetri.f"
	    i__1 = 2, i__2 = ilaenv_(&c__2, "CGETRI", " ", n, &c_n1, &c_n1, &
		    c_n1, (ftnlen)6, (ftnlen)1);
#line 193 "cgetri.f"
	    nbmin = max(i__1,i__2);
#line 194 "cgetri.f"
	}
#line 195 "cgetri.f"
    } else {
#line 196 "cgetri.f"
	iws = *n;
#line 197 "cgetri.f"
    }

/*     Solve the equation inv(A)*L = inv(U) for inv(A). */

#line 201 "cgetri.f"
    if (nb < nbmin || nb >= *n) {

/*        Use unblocked code. */

#line 205 "cgetri.f"
	for (j = *n; j >= 1; --j) {

/*           Copy current column of L to WORK and replace with zeros. */

#line 209 "cgetri.f"
	    i__1 = *n;
#line 209 "cgetri.f"
	    for (i__ = j + 1; i__ <= i__1; ++i__) {
#line 210 "cgetri.f"
		i__2 = i__;
#line 210 "cgetri.f"
		i__3 = i__ + j * a_dim1;
#line 210 "cgetri.f"
		work[i__2].r = a[i__3].r, work[i__2].i = a[i__3].i;
#line 211 "cgetri.f"
		i__2 = i__ + j * a_dim1;
#line 211 "cgetri.f"
		a[i__2].r = 0., a[i__2].i = 0.;
#line 212 "cgetri.f"
/* L10: */
#line 212 "cgetri.f"
	    }

/*           Compute current column of inv(A). */

#line 216 "cgetri.f"
	    if (j < *n) {
#line 216 "cgetri.f"
		i__1 = *n - j;
#line 216 "cgetri.f"
		z__1.r = -1., z__1.i = -0.;
#line 216 "cgetri.f"
		cgemv_("No transpose", n, &i__1, &z__1, &a[(j + 1) * a_dim1 + 
			1], lda, &work[j + 1], &c__1, &c_b2, &a[j * a_dim1 + 
			1], &c__1, (ftnlen)12);
#line 216 "cgetri.f"
	    }
#line 219 "cgetri.f"
/* L20: */
#line 219 "cgetri.f"
	}
#line 220 "cgetri.f"
    } else {

/*        Use blocked code. */

#line 224 "cgetri.f"
	nn = (*n - 1) / nb * nb + 1;
#line 225 "cgetri.f"
	i__1 = -nb;
#line 225 "cgetri.f"
	for (j = nn; i__1 < 0 ? j >= 1 : j <= 1; j += i__1) {
/* Computing MIN */
#line 226 "cgetri.f"
	    i__2 = nb, i__3 = *n - j + 1;
#line 226 "cgetri.f"
	    jb = min(i__2,i__3);

/*           Copy current block column of L to WORK and replace with */
/*           zeros. */

#line 231 "cgetri.f"
	    i__2 = j + jb - 1;
#line 231 "cgetri.f"
	    for (jj = j; jj <= i__2; ++jj) {
#line 232 "cgetri.f"
		i__3 = *n;
#line 232 "cgetri.f"
		for (i__ = jj + 1; i__ <= i__3; ++i__) {
#line 233 "cgetri.f"
		    i__4 = i__ + (jj - j) * ldwork;
#line 233 "cgetri.f"
		    i__5 = i__ + jj * a_dim1;
#line 233 "cgetri.f"
		    work[i__4].r = a[i__5].r, work[i__4].i = a[i__5].i;
#line 234 "cgetri.f"
		    i__4 = i__ + jj * a_dim1;
#line 234 "cgetri.f"
		    a[i__4].r = 0., a[i__4].i = 0.;
#line 235 "cgetri.f"
/* L30: */
#line 235 "cgetri.f"
		}
#line 236 "cgetri.f"
/* L40: */
#line 236 "cgetri.f"
	    }

/*           Compute current block column of inv(A). */

#line 240 "cgetri.f"
	    if (j + jb <= *n) {
#line 240 "cgetri.f"
		i__2 = *n - j - jb + 1;
#line 240 "cgetri.f"
		z__1.r = -1., z__1.i = -0.;
#line 240 "cgetri.f"
		cgemm_("No transpose", "No transpose", n, &jb, &i__2, &z__1, &
			a[(j + jb) * a_dim1 + 1], lda, &work[j + jb], &ldwork,
			 &c_b2, &a[j * a_dim1 + 1], lda, (ftnlen)12, (ftnlen)
			12);
#line 240 "cgetri.f"
	    }
#line 244 "cgetri.f"
	    ctrsm_("Right", "Lower", "No transpose", "Unit", n, &jb, &c_b2, &
		    work[j], &ldwork, &a[j * a_dim1 + 1], lda, (ftnlen)5, (
		    ftnlen)5, (ftnlen)12, (ftnlen)4);
#line 246 "cgetri.f"
/* L50: */
#line 246 "cgetri.f"
	}
#line 247 "cgetri.f"
    }

/*     Apply column interchanges. */

#line 251 "cgetri.f"
    for (j = *n - 1; j >= 1; --j) {
#line 252 "cgetri.f"
	jp = ipiv[j];
#line 253 "cgetri.f"
	if (jp != j) {
#line 253 "cgetri.f"
	    cswap_(n, &a[j * a_dim1 + 1], &c__1, &a[jp * a_dim1 + 1], &c__1);
#line 253 "cgetri.f"
	}
#line 255 "cgetri.f"
/* L60: */
#line 255 "cgetri.f"
    }

#line 257 "cgetri.f"
    work[1].r = (doublereal) iws, work[1].i = 0.;
#line 258 "cgetri.f"
    return 0;

/*     End of CGETRI */

} /* cgetri_ */

