#line 1 "zgeqr.f"
/* zgeqr.f -- translated by f2c (version 20100827).
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

#line 1 "zgeqr.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;


/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGEQR( M, N, A, LDA, T, TSIZE, WORK, LWORK, */
/*                         INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER           INFO, LDA, M, N, TSIZE, LWORK */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16        A( LDA, * ), T( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > ZGEQR computes a QR factorization of an M-by-N matrix A. */
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
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
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
/* >          T is COMPLEX*16 array, dimension (MAX(5,TSIZE)) */
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
/* >          (workspace) COMPLEX*16 array, dimension (MAX(1,LWORK)) */
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
/* >                           ZLATSQR or ZGEQRT */
/* > */
/* >  Depending on the matrix dimensions M and N, and row and column */
/* >  block sizes MB and NB returned by ILAENV, ZGEQR will use either */
/* >  ZLATSQR (if the matrix is tall-and-skinny) or ZGEQRT to compute */
/* >  the QR factorization. */
/* > */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zgeqr_(integer *m, integer *n, doublecomplex *a, integer 
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
    static logical lminws;
    extern /* Subroutine */ int zgeqrt_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static logical lquery;
    static integer mintsz;
    extern /* Subroutine */ int zlatsqr_(integer *, integer *, integer *, 
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

#line 200 "zgeqr.f"
    /* Parameter adjustments */
#line 200 "zgeqr.f"
    a_dim1 = *lda;
#line 200 "zgeqr.f"
    a_offset = 1 + a_dim1;
#line 200 "zgeqr.f"
    a -= a_offset;
#line 200 "zgeqr.f"
    --t;
#line 200 "zgeqr.f"
    --work;
#line 200 "zgeqr.f"

#line 200 "zgeqr.f"
    /* Function Body */
#line 200 "zgeqr.f"
    *info = 0;

#line 202 "zgeqr.f"
    lquery = *tsize == -1 || *tsize == -2 || *lwork == -1 || *lwork == -2;

#line 205 "zgeqr.f"
    mint = FALSE_;
#line 206 "zgeqr.f"
    minw = FALSE_;
#line 207 "zgeqr.f"
    if (*tsize == -2 || *lwork == -2) {
#line 208 "zgeqr.f"
	if (*tsize != -1) {
#line 208 "zgeqr.f"
	    mint = TRUE_;
#line 208 "zgeqr.f"
	}
#line 209 "zgeqr.f"
	if (*lwork != -1) {
#line 209 "zgeqr.f"
	    minw = TRUE_;
#line 209 "zgeqr.f"
	}
#line 210 "zgeqr.f"
    }

/*     Determine the block size */

#line 214 "zgeqr.f"
    if (min(*m,*n) > 0) {
#line 215 "zgeqr.f"
	mb = ilaenv_(&c__1, "ZGEQR ", " ", m, n, &c__1, &c_n1, (ftnlen)6, (
		ftnlen)1);
#line 216 "zgeqr.f"
	nb = ilaenv_(&c__1, "ZGEQR ", " ", m, n, &c__2, &c_n1, (ftnlen)6, (
		ftnlen)1);
#line 217 "zgeqr.f"
    } else {
#line 218 "zgeqr.f"
	mb = *m;
#line 219 "zgeqr.f"
	nb = 1;
#line 220 "zgeqr.f"
    }
#line 221 "zgeqr.f"
    if (mb > *m || mb <= *n) {
#line 221 "zgeqr.f"
	mb = *m;
#line 221 "zgeqr.f"
    }
#line 222 "zgeqr.f"
    if (nb > min(*m,*n) || nb < 1) {
#line 222 "zgeqr.f"
	nb = 1;
#line 222 "zgeqr.f"
    }
#line 223 "zgeqr.f"
    mintsz = *n + 5;
#line 224 "zgeqr.f"
    if (mb > *n && *m > *n) {
#line 225 "zgeqr.f"
	if ((*m - *n) % (mb - *n) == 0) {
#line 226 "zgeqr.f"
	    nblcks = (*m - *n) / (mb - *n);
#line 227 "zgeqr.f"
	} else {
#line 228 "zgeqr.f"
	    nblcks = (*m - *n) / (mb - *n) + 1;
#line 229 "zgeqr.f"
	}
#line 230 "zgeqr.f"
    } else {
#line 231 "zgeqr.f"
	nblcks = 1;
#line 232 "zgeqr.f"
    }

/*     Determine if the workspace size satisfies minimal size */

#line 236 "zgeqr.f"
    lminws = FALSE_;
/* Computing MAX */
#line 237 "zgeqr.f"
    i__1 = 1, i__2 = nb * *n * nblcks + 5;
#line 237 "zgeqr.f"
    if ((*tsize < max(i__1,i__2) || *lwork < nb * *n) && *lwork >= *n && *
	    tsize >= mintsz && ! lquery) {
/* Computing MAX */
#line 240 "zgeqr.f"
	i__1 = 1, i__2 = nb * *n * nblcks + 5;
#line 240 "zgeqr.f"
	if (*tsize < max(i__1,i__2)) {
#line 241 "zgeqr.f"
	    lminws = TRUE_;
#line 242 "zgeqr.f"
	    nb = 1;
#line 243 "zgeqr.f"
	    mb = *m;
#line 244 "zgeqr.f"
	}
#line 245 "zgeqr.f"
	if (*lwork < nb * *n) {
#line 246 "zgeqr.f"
	    lminws = TRUE_;
#line 247 "zgeqr.f"
	    nb = 1;
#line 248 "zgeqr.f"
	}
#line 249 "zgeqr.f"
    }

#line 251 "zgeqr.f"
    if (*m < 0) {
#line 252 "zgeqr.f"
	*info = -1;
#line 253 "zgeqr.f"
    } else if (*n < 0) {
#line 254 "zgeqr.f"
	*info = -2;
#line 255 "zgeqr.f"
    } else if (*lda < max(1,*m)) {
#line 256 "zgeqr.f"
	*info = -4;
#line 257 "zgeqr.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 257 "zgeqr.f"
	i__1 = 1, i__2 = nb * *n * nblcks + 5;
#line 257 "zgeqr.f"
	if (*tsize < max(i__1,i__2) && ! lquery && ! lminws) {
#line 259 "zgeqr.f"
	    *info = -6;
#line 260 "zgeqr.f"
	} else /* if(complicated condition) */ {
/* Computing MAX */
#line 260 "zgeqr.f"
	    i__1 = 1, i__2 = *n * nb;
#line 260 "zgeqr.f"
	    if (*lwork < max(i__1,i__2) && ! lquery && ! lminws) {
#line 262 "zgeqr.f"
		*info = -8;
#line 263 "zgeqr.f"
	    }
#line 263 "zgeqr.f"
	}
#line 263 "zgeqr.f"
    }

#line 265 "zgeqr.f"
    if (*info == 0) {
#line 266 "zgeqr.f"
	if (mint) {
#line 267 "zgeqr.f"
	    t[1].r = (doublereal) mintsz, t[1].i = 0.;
#line 268 "zgeqr.f"
	} else {
#line 269 "zgeqr.f"
	    i__1 = nb * *n * nblcks + 5;
#line 269 "zgeqr.f"
	    t[1].r = (doublereal) i__1, t[1].i = 0.;
#line 270 "zgeqr.f"
	}
#line 271 "zgeqr.f"
	t[2].r = (doublereal) mb, t[2].i = 0.;
#line 272 "zgeqr.f"
	t[3].r = (doublereal) nb, t[3].i = 0.;
#line 273 "zgeqr.f"
	if (minw) {
#line 274 "zgeqr.f"
	    i__1 = max(1,*n);
#line 274 "zgeqr.f"
	    work[1].r = (doublereal) i__1, work[1].i = 0.;
#line 275 "zgeqr.f"
	} else {
/* Computing MAX */
#line 276 "zgeqr.f"
	    i__2 = 1, i__3 = nb * *n;
#line 276 "zgeqr.f"
	    i__1 = max(i__2,i__3);
#line 276 "zgeqr.f"
	    work[1].r = (doublereal) i__1, work[1].i = 0.;
#line 277 "zgeqr.f"
	}
#line 278 "zgeqr.f"
    }
#line 279 "zgeqr.f"
    if (*info != 0) {
#line 280 "zgeqr.f"
	i__1 = -(*info);
#line 280 "zgeqr.f"
	xerbla_("ZGEQR", &i__1, (ftnlen)5);
#line 281 "zgeqr.f"
	return 0;
#line 282 "zgeqr.f"
    } else if (lquery) {
#line 283 "zgeqr.f"
	return 0;
#line 284 "zgeqr.f"
    }

/*     Quick return if possible */

#line 288 "zgeqr.f"
    if (min(*m,*n) == 0) {
#line 289 "zgeqr.f"
	return 0;
#line 290 "zgeqr.f"
    }

/*     The QR Decomposition */

#line 294 "zgeqr.f"
    if (*m <= *n || mb <= *n || mb >= *m) {
#line 295 "zgeqr.f"
	zgeqrt_(m, n, &nb, &a[a_offset], lda, &t[6], &nb, &work[1], info);
#line 296 "zgeqr.f"
    } else {
#line 297 "zgeqr.f"
	zlatsqr_(m, n, &mb, &nb, &a[a_offset], lda, &t[6], &nb, &work[1], 
		lwork, info);
#line 299 "zgeqr.f"
    }

/* Computing MAX */
#line 301 "zgeqr.f"
    i__2 = 1, i__3 = nb * *n;
#line 301 "zgeqr.f"
    i__1 = max(i__2,i__3);
#line 301 "zgeqr.f"
    work[1].r = (doublereal) i__1, work[1].i = 0.;

#line 303 "zgeqr.f"
    return 0;

/*     End of ZGEQR */

} /* zgeqr_ */

