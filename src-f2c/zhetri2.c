#line 1 "zhetri2.f"
/* zhetri2.f -- translated by f2c (version 20100827).
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

#line 1 "zhetri2.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;

/* > \brief \b ZHETRI2 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZHETRI2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetri2
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetri2
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetri2
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHETRI2( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX*16         A( LDA, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHETRI2 computes the inverse of a COMPLEX*16 hermitian indefinite matrix */
/* > A using the factorization A = U*D*U**T or A = L*D*L**T computed by */
/* > ZHETRF. ZHETRI2 set the LEADING DIMENSION of the workspace */
/* > before calling ZHETRI2X that actually computes the inverse. */
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
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the NB diagonal matrix D and the multipliers */
/* >          used to obtain the factor U or L as computed by ZHETRF. */
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
/* >          as determined by ZHETRF. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension (N+NB+1)*(NB+3) */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. */
/* >          WORK is size >= (N+NB+1)*(NB+3) */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >           calculates: */
/* >              - the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, */
/* >              - and no error message related to LWORK is issued by XERBLA. */
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

/* > \date December 2016 */

/* > \ingroup complex16HEcomputational */

/*  ===================================================================== */
/* Subroutine */ int zhetri2_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, doublecomplex *work, integer *lwork, 
	integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;

    /* Local variables */
    extern /* Subroutine */ int zhetri2x_(char *, integer *, doublecomplex *, 
	    integer *, integer *, doublecomplex *, integer *, integer *, 
	    ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nbmax;
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int zhetri_(char *, integer *, doublecomplex *, 
	    integer *, integer *, doublecomplex *, integer *, ftnlen);
    static logical lquery;
    static integer minsize;


/*  -- LAPACK computational routine (version 3.7.0) -- */
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
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 162 "zhetri2.f"
    /* Parameter adjustments */
#line 162 "zhetri2.f"
    a_dim1 = *lda;
#line 162 "zhetri2.f"
    a_offset = 1 + a_dim1;
#line 162 "zhetri2.f"
    a -= a_offset;
#line 162 "zhetri2.f"
    --ipiv;
#line 162 "zhetri2.f"
    --work;
#line 162 "zhetri2.f"

#line 162 "zhetri2.f"
    /* Function Body */
#line 162 "zhetri2.f"
    *info = 0;
#line 163 "zhetri2.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 164 "zhetri2.f"
    lquery = *lwork == -1;
/*     Get blocksize */
#line 166 "zhetri2.f"
    nbmax = ilaenv_(&c__1, "ZHETRF", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6, 
	    (ftnlen)1);
#line 167 "zhetri2.f"
    if (nbmax >= *n) {
#line 168 "zhetri2.f"
	minsize = *n;
#line 169 "zhetri2.f"
    } else {
#line 170 "zhetri2.f"
	minsize = (*n + nbmax + 1) * (nbmax + 3);
#line 171 "zhetri2.f"
    }

#line 173 "zhetri2.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 174 "zhetri2.f"
	*info = -1;
#line 175 "zhetri2.f"
    } else if (*n < 0) {
#line 176 "zhetri2.f"
	*info = -2;
#line 177 "zhetri2.f"
    } else if (*lda < max(1,*n)) {
#line 178 "zhetri2.f"
	*info = -4;
#line 179 "zhetri2.f"
    } else if (*lwork < minsize && ! lquery) {
#line 180 "zhetri2.f"
	*info = -7;
#line 181 "zhetri2.f"
    }

/*     Quick return if possible */


#line 186 "zhetri2.f"
    if (*info != 0) {
#line 187 "zhetri2.f"
	i__1 = -(*info);
#line 187 "zhetri2.f"
	xerbla_("ZHETRI2", &i__1, (ftnlen)7);
#line 188 "zhetri2.f"
	return 0;
#line 189 "zhetri2.f"
    } else if (lquery) {
#line 190 "zhetri2.f"
	work[1].r = (doublereal) minsize, work[1].i = 0.;
#line 191 "zhetri2.f"
	return 0;
#line 192 "zhetri2.f"
    }
#line 193 "zhetri2.f"
    if (*n == 0) {
#line 193 "zhetri2.f"
	return 0;
#line 193 "zhetri2.f"
    }
#line 196 "zhetri2.f"
    if (nbmax >= *n) {
#line 197 "zhetri2.f"
	zhetri_(uplo, n, &a[a_offset], lda, &ipiv[1], &work[1], info, (ftnlen)
		1);
#line 198 "zhetri2.f"
    } else {
#line 199 "zhetri2.f"
	zhetri2x_(uplo, n, &a[a_offset], lda, &ipiv[1], &work[1], &nbmax, 
		info, (ftnlen)1);
#line 200 "zhetri2.f"
    }
#line 201 "zhetri2.f"
    return 0;

/*     End of ZHETRI2 */

} /* zhetri2_ */

