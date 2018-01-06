#line 1 "dgeqr.f"
/* dgeqr.f -- translated by f2c (version 20100827).
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

#line 1 "dgeqr.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;


/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGEQR( M, N, A, LDA, T, TSIZE, WORK, LWORK, */
/*                         INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER           INFO, LDA, M, N, TSIZE, LWORK */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION  A( LDA, * ), T( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > DGEQR computes a QR factorization of an M-by-N matrix A. */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
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
/* >          T is DOUBLE PRECISION array, dimension (MAX(5,TSIZE)) */
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
/* >          (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
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
/* >                           DLATSQR or DGEQRT */
/* > */
/* >  Depending on the matrix dimensions M and N, and row and column */
/* >  block sizes MB and NB returned by ILAENV, DGEQR will use either */
/* >  DLATSQR (if the matrix is tall-and-skinny) or DGEQRT to compute */
/* >  the QR factorization. */
/* > */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dgeqr_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *t, integer *tsize, doublereal *work, integer *lwork, 
	integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer mb, nb;
    static logical mint, minw;
    static integer nblcks;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgeqrt_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static logical lminws, lquery;
    static integer mintsz;
    extern /* Subroutine */ int dlatsqr_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *);


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

#line 200 "dgeqr.f"
    /* Parameter adjustments */
#line 200 "dgeqr.f"
    a_dim1 = *lda;
#line 200 "dgeqr.f"
    a_offset = 1 + a_dim1;
#line 200 "dgeqr.f"
    a -= a_offset;
#line 200 "dgeqr.f"
    --t;
#line 200 "dgeqr.f"
    --work;
#line 200 "dgeqr.f"

#line 200 "dgeqr.f"
    /* Function Body */
#line 200 "dgeqr.f"
    *info = 0;

#line 202 "dgeqr.f"
    lquery = *tsize == -1 || *tsize == -2 || *lwork == -1 || *lwork == -2;

#line 205 "dgeqr.f"
    mint = FALSE_;
#line 206 "dgeqr.f"
    minw = FALSE_;
#line 207 "dgeqr.f"
    if (*tsize == -2 || *lwork == -2) {
#line 208 "dgeqr.f"
	if (*tsize != -1) {
#line 208 "dgeqr.f"
	    mint = TRUE_;
#line 208 "dgeqr.f"
	}
#line 209 "dgeqr.f"
	if (*lwork != -1) {
#line 209 "dgeqr.f"
	    minw = TRUE_;
#line 209 "dgeqr.f"
	}
#line 210 "dgeqr.f"
    }

/*     Determine the block size */

#line 214 "dgeqr.f"
    if (min(*m,*n) > 0) {
#line 215 "dgeqr.f"
	mb = ilaenv_(&c__1, "DGEQR ", " ", m, n, &c__1, &c_n1, (ftnlen)6, (
		ftnlen)1);
#line 216 "dgeqr.f"
	nb = ilaenv_(&c__1, "DGEQR ", " ", m, n, &c__2, &c_n1, (ftnlen)6, (
		ftnlen)1);
#line 217 "dgeqr.f"
    } else {
#line 218 "dgeqr.f"
	mb = *m;
#line 219 "dgeqr.f"
	nb = 1;
#line 220 "dgeqr.f"
    }
#line 221 "dgeqr.f"
    if (mb > *m || mb <= *n) {
#line 221 "dgeqr.f"
	mb = *m;
#line 221 "dgeqr.f"
    }
#line 222 "dgeqr.f"
    if (nb > min(*m,*n) || nb < 1) {
#line 222 "dgeqr.f"
	nb = 1;
#line 222 "dgeqr.f"
    }
#line 223 "dgeqr.f"
    mintsz = *n + 5;
#line 224 "dgeqr.f"
    if (mb > *n && *m > *n) {
#line 225 "dgeqr.f"
	if ((*m - *n) % (mb - *n) == 0) {
#line 226 "dgeqr.f"
	    nblcks = (*m - *n) / (mb - *n);
#line 227 "dgeqr.f"
	} else {
#line 228 "dgeqr.f"
	    nblcks = (*m - *n) / (mb - *n) + 1;
#line 229 "dgeqr.f"
	}
#line 230 "dgeqr.f"
    } else {
#line 231 "dgeqr.f"
	nblcks = 1;
#line 232 "dgeqr.f"
    }

/*     Determine if the workspace size satisfies minimal size */

#line 236 "dgeqr.f"
    lminws = FALSE_;
/* Computing MAX */
#line 237 "dgeqr.f"
    i__1 = 1, i__2 = nb * *n * nblcks + 5;
#line 237 "dgeqr.f"
    if ((*tsize < max(i__1,i__2) || *lwork < nb * *n) && *lwork >= *n && *
	    tsize >= mintsz && ! lquery) {
/* Computing MAX */
#line 240 "dgeqr.f"
	i__1 = 1, i__2 = nb * *n * nblcks + 5;
#line 240 "dgeqr.f"
	if (*tsize < max(i__1,i__2)) {
#line 241 "dgeqr.f"
	    lminws = TRUE_;
#line 242 "dgeqr.f"
	    nb = 1;
#line 243 "dgeqr.f"
	    mb = *m;
#line 244 "dgeqr.f"
	}
#line 245 "dgeqr.f"
	if (*lwork < nb * *n) {
#line 246 "dgeqr.f"
	    lminws = TRUE_;
#line 247 "dgeqr.f"
	    nb = 1;
#line 248 "dgeqr.f"
	}
#line 249 "dgeqr.f"
    }

#line 251 "dgeqr.f"
    if (*m < 0) {
#line 252 "dgeqr.f"
	*info = -1;
#line 253 "dgeqr.f"
    } else if (*n < 0) {
#line 254 "dgeqr.f"
	*info = -2;
#line 255 "dgeqr.f"
    } else if (*lda < max(1,*m)) {
#line 256 "dgeqr.f"
	*info = -4;
#line 257 "dgeqr.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 257 "dgeqr.f"
	i__1 = 1, i__2 = nb * *n * nblcks + 5;
#line 257 "dgeqr.f"
	if (*tsize < max(i__1,i__2) && ! lquery && ! lminws) {
#line 259 "dgeqr.f"
	    *info = -6;
#line 260 "dgeqr.f"
	} else /* if(complicated condition) */ {
/* Computing MAX */
#line 260 "dgeqr.f"
	    i__1 = 1, i__2 = *n * nb;
#line 260 "dgeqr.f"
	    if (*lwork < max(i__1,i__2) && ! lquery && ! lminws) {
#line 262 "dgeqr.f"
		*info = -8;
#line 263 "dgeqr.f"
	    }
#line 263 "dgeqr.f"
	}
#line 263 "dgeqr.f"
    }

#line 265 "dgeqr.f"
    if (*info == 0) {
#line 266 "dgeqr.f"
	if (mint) {
#line 267 "dgeqr.f"
	    t[1] = (doublereal) mintsz;
#line 268 "dgeqr.f"
	} else {
#line 269 "dgeqr.f"
	    t[1] = (doublereal) (nb * *n * nblcks + 5);
#line 270 "dgeqr.f"
	}
#line 271 "dgeqr.f"
	t[2] = (doublereal) mb;
#line 272 "dgeqr.f"
	t[3] = (doublereal) nb;
#line 273 "dgeqr.f"
	if (minw) {
#line 274 "dgeqr.f"
	    work[1] = (doublereal) max(1,*n);
#line 275 "dgeqr.f"
	} else {
/* Computing MAX */
#line 276 "dgeqr.f"
	    i__1 = 1, i__2 = nb * *n;
#line 276 "dgeqr.f"
	    work[1] = (doublereal) max(i__1,i__2);
#line 277 "dgeqr.f"
	}
#line 278 "dgeqr.f"
    }
#line 279 "dgeqr.f"
    if (*info != 0) {
#line 280 "dgeqr.f"
	i__1 = -(*info);
#line 280 "dgeqr.f"
	xerbla_("DGEQR", &i__1, (ftnlen)5);
#line 281 "dgeqr.f"
	return 0;
#line 282 "dgeqr.f"
    } else if (lquery) {
#line 283 "dgeqr.f"
	return 0;
#line 284 "dgeqr.f"
    }

/*     Quick return if possible */

#line 288 "dgeqr.f"
    if (min(*m,*n) == 0) {
#line 289 "dgeqr.f"
	return 0;
#line 290 "dgeqr.f"
    }

/*     The QR Decomposition */

#line 294 "dgeqr.f"
    if (*m <= *n || mb <= *n || mb >= *m) {
#line 295 "dgeqr.f"
	dgeqrt_(m, n, &nb, &a[a_offset], lda, &t[6], &nb, &work[1], info);
#line 296 "dgeqr.f"
    } else {
#line 297 "dgeqr.f"
	dlatsqr_(m, n, &mb, &nb, &a[a_offset], lda, &t[6], &nb, &work[1], 
		lwork, info);
#line 299 "dgeqr.f"
    }

/* Computing MAX */
#line 301 "dgeqr.f"
    i__1 = 1, i__2 = nb * *n;
#line 301 "dgeqr.f"
    work[1] = (doublereal) max(i__1,i__2);

#line 303 "dgeqr.f"
    return 0;

/*     End of DGEQR */

} /* dgeqr_ */

