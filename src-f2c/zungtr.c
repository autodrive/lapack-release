#line 1 "zungtr.f"
/* zungtr.f -- translated by f2c (version 20100827).
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

#line 1 "zungtr.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;

/* > \brief \b ZUNGTR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZUNGTR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zungtr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zungtr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zungtr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZUNGTR( UPLO, N, A, LDA, TAU, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZUNGTR generates a complex unitary matrix Q which is defined as the */
/* > product of n-1 elementary reflectors of order N, as returned by */
/* > ZHETRD: */
/* > */
/* > if UPLO = 'U', Q = H(n-1) . . . H(2) H(1), */
/* > */
/* > if UPLO = 'L', Q = H(1) H(2) . . . H(n-1). */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U': Upper triangle of A contains elementary reflectors */
/* >                 from ZHETRD; */
/* >          = 'L': Lower triangle of A contains elementary reflectors */
/* >                 from ZHETRD. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix Q. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the vectors which define the elementary reflectors, */
/* >          as returned by ZHETRD. */
/* >          On exit, the N-by-N unitary matrix Q. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A. LDA >= N. */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* >          TAU is COMPLEX*16 array, dimension (N-1) */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i), as returned by ZHETRD. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. LWORK >= N-1. */
/* >          For optimum performance LWORK >= (N-1)*NB, where NB is */
/* >          the optimal blocksize. */
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
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complex16OTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int zungtr_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, doublecomplex *work, integer *lwork,
	 integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, j, nb;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer lwkopt;
    static logical lquery;
    extern /* Subroutine */ int zungql_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *), zungqr_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *);


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

/*     Test the input arguments */

#line 165 "zungtr.f"
    /* Parameter adjustments */
#line 165 "zungtr.f"
    a_dim1 = *lda;
#line 165 "zungtr.f"
    a_offset = 1 + a_dim1;
#line 165 "zungtr.f"
    a -= a_offset;
#line 165 "zungtr.f"
    --tau;
#line 165 "zungtr.f"
    --work;
#line 165 "zungtr.f"

#line 165 "zungtr.f"
    /* Function Body */
#line 165 "zungtr.f"
    *info = 0;
#line 166 "zungtr.f"
    lquery = *lwork == -1;
#line 167 "zungtr.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 168 "zungtr.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 169 "zungtr.f"
	*info = -1;
#line 170 "zungtr.f"
    } else if (*n < 0) {
#line 171 "zungtr.f"
	*info = -2;
#line 172 "zungtr.f"
    } else if (*lda < max(1,*n)) {
#line 173 "zungtr.f"
	*info = -4;
#line 174 "zungtr.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 174 "zungtr.f"
	i__1 = 1, i__2 = *n - 1;
#line 174 "zungtr.f"
	if (*lwork < max(i__1,i__2) && ! lquery) {
#line 175 "zungtr.f"
	    *info = -7;
#line 176 "zungtr.f"
	}
#line 176 "zungtr.f"
    }

#line 178 "zungtr.f"
    if (*info == 0) {
#line 179 "zungtr.f"
	if (upper) {
#line 180 "zungtr.f"
	    i__1 = *n - 1;
#line 180 "zungtr.f"
	    i__2 = *n - 1;
#line 180 "zungtr.f"
	    i__3 = *n - 1;
#line 180 "zungtr.f"
	    nb = ilaenv_(&c__1, "ZUNGQL", " ", &i__1, &i__2, &i__3, &c_n1, (
		    ftnlen)6, (ftnlen)1);
#line 181 "zungtr.f"
	} else {
#line 182 "zungtr.f"
	    i__1 = *n - 1;
#line 182 "zungtr.f"
	    i__2 = *n - 1;
#line 182 "zungtr.f"
	    i__3 = *n - 1;
#line 182 "zungtr.f"
	    nb = ilaenv_(&c__1, "ZUNGQR", " ", &i__1, &i__2, &i__3, &c_n1, (
		    ftnlen)6, (ftnlen)1);
#line 183 "zungtr.f"
	}
/* Computing MAX */
#line 184 "zungtr.f"
	i__1 = 1, i__2 = *n - 1;
#line 184 "zungtr.f"
	lwkopt = max(i__1,i__2) * nb;
#line 185 "zungtr.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 186 "zungtr.f"
    }

#line 188 "zungtr.f"
    if (*info != 0) {
#line 189 "zungtr.f"
	i__1 = -(*info);
#line 189 "zungtr.f"
	xerbla_("ZUNGTR", &i__1, (ftnlen)6);
#line 190 "zungtr.f"
	return 0;
#line 191 "zungtr.f"
    } else if (lquery) {
#line 192 "zungtr.f"
	return 0;
#line 193 "zungtr.f"
    }

/*     Quick return if possible */

#line 197 "zungtr.f"
    if (*n == 0) {
#line 198 "zungtr.f"
	work[1].r = 1., work[1].i = 0.;
#line 199 "zungtr.f"
	return 0;
#line 200 "zungtr.f"
    }

#line 202 "zungtr.f"
    if (upper) {

/*        Q was determined by a call to ZHETRD with UPLO = 'U' */

/*        Shift the vectors which define the elementary reflectors one */
/*        column to the left, and set the last row and column of Q to */
/*        those of the unit matrix */

#line 210 "zungtr.f"
	i__1 = *n - 1;
#line 210 "zungtr.f"
	for (j = 1; j <= i__1; ++j) {
#line 211 "zungtr.f"
	    i__2 = j - 1;
#line 211 "zungtr.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 212 "zungtr.f"
		i__3 = i__ + j * a_dim1;
#line 212 "zungtr.f"
		i__4 = i__ + (j + 1) * a_dim1;
#line 212 "zungtr.f"
		a[i__3].r = a[i__4].r, a[i__3].i = a[i__4].i;
#line 213 "zungtr.f"
/* L10: */
#line 213 "zungtr.f"
	    }
#line 214 "zungtr.f"
	    i__2 = *n + j * a_dim1;
#line 214 "zungtr.f"
	    a[i__2].r = 0., a[i__2].i = 0.;
#line 215 "zungtr.f"
/* L20: */
#line 215 "zungtr.f"
	}
#line 216 "zungtr.f"
	i__1 = *n - 1;
#line 216 "zungtr.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 217 "zungtr.f"
	    i__2 = i__ + *n * a_dim1;
#line 217 "zungtr.f"
	    a[i__2].r = 0., a[i__2].i = 0.;
#line 218 "zungtr.f"
/* L30: */
#line 218 "zungtr.f"
	}
#line 219 "zungtr.f"
	i__1 = *n + *n * a_dim1;
#line 219 "zungtr.f"
	a[i__1].r = 1., a[i__1].i = 0.;

/*        Generate Q(1:n-1,1:n-1) */

#line 223 "zungtr.f"
	i__1 = *n - 1;
#line 223 "zungtr.f"
	i__2 = *n - 1;
#line 223 "zungtr.f"
	i__3 = *n - 1;
#line 223 "zungtr.f"
	zungql_(&i__1, &i__2, &i__3, &a[a_offset], lda, &tau[1], &work[1], 
		lwork, &iinfo);

#line 225 "zungtr.f"
    } else {

/*        Q was determined by a call to ZHETRD with UPLO = 'L'. */

/*        Shift the vectors which define the elementary reflectors one */
/*        column to the right, and set the first row and column of Q to */
/*        those of the unit matrix */

#line 233 "zungtr.f"
	for (j = *n; j >= 2; --j) {
#line 234 "zungtr.f"
	    i__1 = j * a_dim1 + 1;
#line 234 "zungtr.f"
	    a[i__1].r = 0., a[i__1].i = 0.;
#line 235 "zungtr.f"
	    i__1 = *n;
#line 235 "zungtr.f"
	    for (i__ = j + 1; i__ <= i__1; ++i__) {
#line 236 "zungtr.f"
		i__2 = i__ + j * a_dim1;
#line 236 "zungtr.f"
		i__3 = i__ + (j - 1) * a_dim1;
#line 236 "zungtr.f"
		a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
#line 237 "zungtr.f"
/* L40: */
#line 237 "zungtr.f"
	    }
#line 238 "zungtr.f"
/* L50: */
#line 238 "zungtr.f"
	}
#line 239 "zungtr.f"
	i__1 = a_dim1 + 1;
#line 239 "zungtr.f"
	a[i__1].r = 1., a[i__1].i = 0.;
#line 240 "zungtr.f"
	i__1 = *n;
#line 240 "zungtr.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 241 "zungtr.f"
	    i__2 = i__ + a_dim1;
#line 241 "zungtr.f"
	    a[i__2].r = 0., a[i__2].i = 0.;
#line 242 "zungtr.f"
/* L60: */
#line 242 "zungtr.f"
	}
#line 243 "zungtr.f"
	if (*n > 1) {

/*           Generate Q(2:n,2:n) */

#line 247 "zungtr.f"
	    i__1 = *n - 1;
#line 247 "zungtr.f"
	    i__2 = *n - 1;
#line 247 "zungtr.f"
	    i__3 = *n - 1;
#line 247 "zungtr.f"
	    zungqr_(&i__1, &i__2, &i__3, &a[(a_dim1 << 1) + 2], lda, &tau[1], 
		    &work[1], lwork, &iinfo);
#line 249 "zungtr.f"
	}
#line 250 "zungtr.f"
    }
#line 251 "zungtr.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 252 "zungtr.f"
    return 0;

/*     End of ZUNGTR */

} /* zungtr_ */

