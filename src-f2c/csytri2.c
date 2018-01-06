#line 1 "csytri2.f"
/* csytri2.f -- translated by f2c (version 20100827).
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

#line 1 "csytri2.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;

/* > \brief \b CSYTRI2 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CSYTRI2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csytri2
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csytri2
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csytri2
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CSYTRI2( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
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
/* > CSYTRI2 computes the inverse of a COMPLEX symmetric indefinite matrix */
/* > A using the factorization A = U*D*U**T or A = L*D*L**T computed by */
/* > CSYTRF. CSYTRI2 sets the LEADING DIMENSION of the workspace */
/* > before calling CSYTRI2X that actually computes the inverse. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the details of the factorization are stored */
/* >          as an upper or lower triangular matrix. */
/* >          = 'U':  Upper triangular, form is A = U*D*U**T; */
/* >          = 'L':  Lower triangular, form is A = L*D*L**T. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          On entry, the NB diagonal matrix D and the multipliers */
/* >          used to obtain the factor U or L as computed by CSYTRF. */
/* > */
/* >          On exit, if INFO = 0, the (symmetric) inverse of the original */
/* >          matrix.  If UPLO = 'U', the upper triangular part of the */
/* >          inverse is formed and the part of A below the diagonal is not */
/* >          referenced; if UPLO = 'L' the lower triangular part of the */
/* >          inverse is formed and the part of A above the diagonal is */
/* >          not referenced. */
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
/* >          Details of the interchanges and the NB structure of D */
/* >          as determined by CSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX array, dimension (N+NB+1)*(NB+3) */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. */
/* >          WORK is size >= (N+NB+1)*(NB+3) */
/* >          If LDWORK = -1, then a workspace query is assumed; the routine */
/* >           calculates: */
/* >              - the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, */
/* >              - and no error message related to LDWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0: successful exit */
/* >          < 0: if INFO = -i, the i-th argument had an illegal value */
/* >          > 0: if INFO = i, D(i,i) = 0; the matrix is singular and its */
/* >               inverse could not be computed. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup complexSYcomputational */

/*  ===================================================================== */
/* Subroutine */ int csytri2_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, doublecomplex *work, integer *lwork, 
	integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;

    /* Local variables */
    extern /* Subroutine */ int csytri2x_(char *, integer *, doublecomplex *, 
	    integer *, integer *, doublecomplex *, integer *, integer *, 
	    ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nbmax;
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int csytri_(char *, integer *, doublecomplex *, 
	    integer *, integer *, doublecomplex *, integer *, ftnlen);
    static logical lquery;
    static integer minsize;


/*  -- LAPACK computational routine (version 3.4.2) -- */
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
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 162 "csytri2.f"
    /* Parameter adjustments */
#line 162 "csytri2.f"
    a_dim1 = *lda;
#line 162 "csytri2.f"
    a_offset = 1 + a_dim1;
#line 162 "csytri2.f"
    a -= a_offset;
#line 162 "csytri2.f"
    --ipiv;
#line 162 "csytri2.f"
    --work;
#line 162 "csytri2.f"

#line 162 "csytri2.f"
    /* Function Body */
#line 162 "csytri2.f"
    *info = 0;
#line 163 "csytri2.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 164 "csytri2.f"
    lquery = *lwork == -1;
/*     Get blocksize */
#line 166 "csytri2.f"
    nbmax = ilaenv_(&c__1, "CSYTRF", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6, 
	    (ftnlen)1);
#line 167 "csytri2.f"
    if (nbmax >= *n) {
#line 168 "csytri2.f"
	minsize = *n;
#line 169 "csytri2.f"
    } else {
#line 170 "csytri2.f"
	minsize = (*n + nbmax + 1) * (nbmax + 3);
#line 171 "csytri2.f"
    }

#line 173 "csytri2.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 174 "csytri2.f"
	*info = -1;
#line 175 "csytri2.f"
    } else if (*n < 0) {
#line 176 "csytri2.f"
	*info = -2;
#line 177 "csytri2.f"
    } else if (*lda < max(1,*n)) {
#line 178 "csytri2.f"
	*info = -4;
#line 179 "csytri2.f"
    } else if (*lwork < minsize && ! lquery) {
#line 180 "csytri2.f"
	*info = -7;
#line 181 "csytri2.f"
    }

/*     Quick return if possible */


#line 186 "csytri2.f"
    if (*info != 0) {
#line 187 "csytri2.f"
	i__1 = -(*info);
#line 187 "csytri2.f"
	xerbla_("CSYTRI2", &i__1, (ftnlen)7);
#line 188 "csytri2.f"
	return 0;
#line 189 "csytri2.f"
    } else if (lquery) {
#line 190 "csytri2.f"
	work[1].r = (doublereal) minsize, work[1].i = 0.;
#line 191 "csytri2.f"
	return 0;
#line 192 "csytri2.f"
    }
#line 193 "csytri2.f"
    if (*n == 0) {
#line 193 "csytri2.f"
	return 0;
#line 193 "csytri2.f"
    }
#line 196 "csytri2.f"
    if (nbmax >= *n) {
#line 197 "csytri2.f"
	csytri_(uplo, n, &a[a_offset], lda, &ipiv[1], &work[1], info, (ftnlen)
		1);
#line 198 "csytri2.f"
    } else {
#line 199 "csytri2.f"
	csytri2x_(uplo, n, &a[a_offset], lda, &ipiv[1], &work[1], &nbmax, 
		info, (ftnlen)1);
#line 200 "csytri2.f"
    }
#line 201 "csytri2.f"
    return 0;

/*     End of CSYTRI2 */

} /* csytri2_ */

