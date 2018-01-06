#line 1 "cgeqr.f"
/* cgeqr.f -- translated by f2c (version 20100827).
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

#line 1 "cgeqr.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;


/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGEQR( M, N, A, LDA, T, TSIZE, WORK, LWORK, */
/*                         INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER           INFO, LDA, M, N, TSIZE, LWORK */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX           A( LDA, * ), T( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > CGEQR computes a QR factorization of an M-by-N matrix A. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix A.  M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, the elements on and above the diagonal of the array */
/* >          contain the min(M,N)-by-N upper trapezoidal matrix R */
/* >          (R is upper triangular if M >= N); */
/* >          the elements below the diagonal are used to store part of the */
/* >          data structure to represent Q. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] T */
/* > \verbatim */
/* >          T is COMPLEX array, dimension (MAX(5,TSIZE)) */
/* >          On exit, if INFO = 0, T(1) returns optimal (or either minimal */
/* >          or optimal, if query is assumed) TSIZE. See TSIZE for details. */
/* >          Remaining T contains part of the data structure used to represent Q. */
/* >          If one wants to apply or construct Q, then one needs to keep T */
/* >          (in addition to A) and pass it to further subroutines. */
/* > \endverbatim */
/* > */
/* > \param[in] TSIZE */
/* > \verbatim */
/* >          TSIZE is INTEGER */
/* >          If TSIZE >= 5, the dimension of the array T. */
/* >          If TSIZE = -1 or -2, then a workspace query is assumed. The routine */
/* >          only calculates the sizes of the T and WORK arrays, returns these */
/* >          values as the first entries of the T and WORK arrays, and no error */
/* >          message related to T or WORK is issued by XERBLA. */
/* >          If TSIZE = -1, the routine calculates optimal size of T for the */
/* >          optimum performance and returns this value in T(1). */
/* >          If TSIZE = -2, the routine calculates minimal size of T and */
/* >          returns this value in T(1). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          (workspace) COMPLEX array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) contains optimal (or either minimal */
/* >          or optimal, if query was assumed) LWORK. */
/* >          See LWORK for details. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. */
/* >          If LWORK = -1 or -2, then a workspace query is assumed. The routine */
/* >          only calculates the sizes of the T and WORK arrays, returns these */
/* >          values as the first entries of the T and WORK arrays, and no error */
/* >          message related to T or WORK is issued by XERBLA. */
/* >          If LWORK = -1, the routine calculates optimal size of WORK for the */
/* >          optimal performance and returns this value in WORK(1). */
/* >          If LWORK = -2, the routine calculates minimal size of WORK and */
/* >          returns this value in WORK(1). */
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

/* > \par Further Details */
/*  ==================== */
/* > */
/* > \verbatim */
/* > */
/* > The goal of the interface is to give maximum freedom to the developers for */
/* > creating any QR factorization algorithm they wish. The triangular */
/* > (trapezoidal) R has to be stored in the upper part of A. The lower part of A */
/* > and the array T can be used to store any relevant information for applying or */
/* > constructing the Q factor. The WORK array can safely be discarded after exit. */
/* > */
/* > Caution: One should not expect the sizes of T and WORK to be the same from one */
/* > LAPACK implementation to the other, or even from one execution to the other. */
/* > A workspace query (for T and WORK) is needed at each execution. However, */
/* > for a given execution, the size of T and WORK are fixed and will not change */
/* > from one query to the next. */
/* > */
/* > \endverbatim */
/* > */
/* > \par Further Details particular to this LAPACK implementation: */
/*  ============================================================== */
/* > */
/* > \verbatim */
/* > */
/* > These details are particular for this LAPACK implementation. Users should not */
/* > take them for granted. These details may change in the future, and are unlikely not */
/* > true for another LAPACK implementation. These details are relevant if one wants */
/* > to try to understand the code. They are not part of the interface. */
/* > */
/* > In this version, */
/* > */
/* >          T(2): row block size (MB) */
/* >          T(3): column block size (NB) */
/* >          T(6:TSIZE): data structure needed for Q, computed by */
/* >                           CLATSQR or CGEQRT */
/* > */
/* >  Depending on the matrix dimensions M and N, and row and column */
/* >  block sizes MB and NB returned by ILAENV, CGEQR will use either */
/* >  CLATSQR (if the matrix is tall-and-skinny) or CGEQRT to compute */
/* >  the QR factorization. */
/* > */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int cgeqr_(integer *m, integer *n, doublecomplex *a, integer 
	*lda, doublecomplex *t, integer *tsize, doublecomplex *work, integer *
	lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer mb, nb;
    static logical mint, minw;
    static integer nblcks;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int cgeqrt_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static logical lminws, lquery;
    static integer mintsz;
    extern /* Subroutine */ int clatsqr_(integer *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublecomplex *, integer *, integer *);


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd. -- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments */

#line 200 "cgeqr.f"
    /* Parameter adjustments */
#line 200 "cgeqr.f"
    a_dim1 = *lda;
#line 200 "cgeqr.f"
    a_offset = 1 + a_dim1;
#line 200 "cgeqr.f"
    a -= a_offset;
#line 200 "cgeqr.f"
    --t;
#line 200 "cgeqr.f"
    --work;
#line 200 "cgeqr.f"

#line 200 "cgeqr.f"
    /* Function Body */
#line 200 "cgeqr.f"
    *info = 0;

#line 202 "cgeqr.f"
    lquery = *tsize == -1 || *tsize == -2 || *lwork == -1 || *lwork == -2;

#line 205 "cgeqr.f"
    mint = FALSE_;
#line 206 "cgeqr.f"
    minw = FALSE_;
#line 207 "cgeqr.f"
    if (*tsize == -2 || *lwork == -2) {
#line 208 "cgeqr.f"
	if (*tsize != -1) {
#line 208 "cgeqr.f"
	    mint = TRUE_;
#line 208 "cgeqr.f"
	}
#line 209 "cgeqr.f"
	if (*lwork != -1) {
#line 209 "cgeqr.f"
	    minw = TRUE_;
#line 209 "cgeqr.f"
	}
#line 210 "cgeqr.f"
    }

/*     Determine the block size */

#line 214 "cgeqr.f"
    if (min(*m,*n) > 0) {
#line 215 "cgeqr.f"
	mb = ilaenv_(&c__1, "CGEQR ", " ", m, n, &c__1, &c_n1, (ftnlen)6, (
		ftnlen)1);
#line 216 "cgeqr.f"
	nb = ilaenv_(&c__1, "CGEQR ", " ", m, n, &c__2, &c_n1, (ftnlen)6, (
		ftnlen)1);
#line 217 "cgeqr.f"
    } else {
#line 218 "cgeqr.f"
	mb = *m;
#line 219 "cgeqr.f"
	nb = 1;
#line 220 "cgeqr.f"
    }
#line 221 "cgeqr.f"
    if (mb > *m || mb <= *n) {
#line 221 "cgeqr.f"
	mb = *m;
#line 221 "cgeqr.f"
    }
#line 222 "cgeqr.f"
    if (nb > min(*m,*n) || nb < 1) {
#line 222 "cgeqr.f"
	nb = 1;
#line 222 "cgeqr.f"
    }
#line 223 "cgeqr.f"
    mintsz = *n + 5;
#line 224 "cgeqr.f"
    if (mb > *n && *m > *n) {
#line 225 "cgeqr.f"
	if ((*m - *n) % (mb - *n) == 0) {
#line 226 "cgeqr.f"
	    nblcks = (*m - *n) / (mb - *n);
#line 227 "cgeqr.f"
	} else {
#line 228 "cgeqr.f"
	    nblcks = (*m - *n) / (mb - *n) + 1;
#line 229 "cgeqr.f"
	}
#line 230 "cgeqr.f"
    } else {
#line 231 "cgeqr.f"
	nblcks = 1;
#line 232 "cgeqr.f"
    }

/*     Determine if the workspace size satisfies minimal size */

#line 236 "cgeqr.f"
    lminws = FALSE_;
/* Computing MAX */
#line 237 "cgeqr.f"
    i__1 = 1, i__2 = nb * *n * nblcks + 5;
#line 237 "cgeqr.f"
    if ((*tsize < max(i__1,i__2) || *lwork < nb * *n) && *lwork >= *n && *
	    tsize >= mintsz && ! lquery) {
/* Computing MAX */
#line 240 "cgeqr.f"
	i__1 = 1, i__2 = nb * *n * nblcks + 5;
#line 240 "cgeqr.f"
	if (*tsize < max(i__1,i__2)) {
#line 241 "cgeqr.f"
	    lminws = TRUE_;
#line 242 "cgeqr.f"
	    nb = 1;
#line 243 "cgeqr.f"
	    mb = *m;
#line 244 "cgeqr.f"
	}
#line 245 "cgeqr.f"
	if (*lwork < nb * *n) {
#line 246 "cgeqr.f"
	    lminws = TRUE_;
#line 247 "cgeqr.f"
	    nb = 1;
#line 248 "cgeqr.f"
	}
#line 249 "cgeqr.f"
    }

#line 251 "cgeqr.f"
    if (*m < 0) {
#line 252 "cgeqr.f"
	*info = -1;
#line 253 "cgeqr.f"
    } else if (*n < 0) {
#line 254 "cgeqr.f"
	*info = -2;
#line 255 "cgeqr.f"
    } else if (*lda < max(1,*m)) {
#line 256 "cgeqr.f"
	*info = -4;
#line 257 "cgeqr.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 257 "cgeqr.f"
	i__1 = 1, i__2 = nb * *n * nblcks + 5;
#line 257 "cgeqr.f"
	if (*tsize < max(i__1,i__2) && ! lquery && ! lminws) {
#line 259 "cgeqr.f"
	    *info = -6;
#line 260 "cgeqr.f"
	} else /* if(complicated condition) */ {
/* Computing MAX */
#line 260 "cgeqr.f"
	    i__1 = 1, i__2 = *n * nb;
#line 260 "cgeqr.f"
	    if (*lwork < max(i__1,i__2) && ! lquery && ! lminws) {
#line 262 "cgeqr.f"
		*info = -8;
#line 263 "cgeqr.f"
	    }
#line 263 "cgeqr.f"
	}
#line 263 "cgeqr.f"
    }

#line 265 "cgeqr.f"
    if (*info == 0) {
#line 266 "cgeqr.f"
	if (mint) {
#line 267 "cgeqr.f"
	    t[1].r = (doublereal) mintsz, t[1].i = 0.;
#line 268 "cgeqr.f"
	} else {
#line 269 "cgeqr.f"
	    i__1 = nb * *n * nblcks + 5;
#line 269 "cgeqr.f"
	    t[1].r = (doublereal) i__1, t[1].i = 0.;
#line 270 "cgeqr.f"
	}
#line 271 "cgeqr.f"
	t[2].r = (doublereal) mb, t[2].i = 0.;
#line 272 "cgeqr.f"
	t[3].r = (doublereal) nb, t[3].i = 0.;
#line 273 "cgeqr.f"
	if (minw) {
#line 274 "cgeqr.f"
	    i__1 = max(1,*n);
#line 274 "cgeqr.f"
	    work[1].r = (doublereal) i__1, work[1].i = 0.;
#line 275 "cgeqr.f"
	} else {
/* Computing MAX */
#line 276 "cgeqr.f"
	    i__2 = 1, i__3 = nb * *n;
#line 276 "cgeqr.f"
	    i__1 = max(i__2,i__3);
#line 276 "cgeqr.f"
	    work[1].r = (doublereal) i__1, work[1].i = 0.;
#line 277 "cgeqr.f"
	}
#line 278 "cgeqr.f"
    }
#line 279 "cgeqr.f"
    if (*info != 0) {
#line 280 "cgeqr.f"
	i__1 = -(*info);
#line 280 "cgeqr.f"
	xerbla_("CGEQR", &i__1, (ftnlen)5);
#line 281 "cgeqr.f"
	return 0;
#line 282 "cgeqr.f"
    } else if (lquery) {
#line 283 "cgeqr.f"
	return 0;
#line 284 "cgeqr.f"
    }

/*     Quick return if possible */

#line 288 "cgeqr.f"
    if (min(*m,*n) == 0) {
#line 289 "cgeqr.f"
	return 0;
#line 290 "cgeqr.f"
    }

/*     The QR Decomposition */

#line 294 "cgeqr.f"
    if (*m <= *n || mb <= *n || mb >= *m) {
#line 295 "cgeqr.f"
	cgeqrt_(m, n, &nb, &a[a_offset], lda, &t[6], &nb, &work[1], info);
#line 296 "cgeqr.f"
    } else {
#line 297 "cgeqr.f"
	clatsqr_(m, n, &mb, &nb, &a[a_offset], lda, &t[6], &nb, &work[1], 
		lwork, info);
#line 299 "cgeqr.f"
    }

/* Computing MAX */
#line 301 "cgeqr.f"
    i__2 = 1, i__3 = nb * *n;
#line 301 "cgeqr.f"
    i__1 = max(i__2,i__3);
#line 301 "cgeqr.f"
    work[1].r = (doublereal) i__1, work[1].i = 0.;

#line 303 "cgeqr.f"
    return 0;

/*     End of CGEQR */

} /* cgeqr_ */

