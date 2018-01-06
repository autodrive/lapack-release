#line 1 "cgelq.f"
/* cgelq.f -- translated by f2c (version 20100827).
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

#line 1 "cgelq.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;


/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGELQ( M, N, A, LDA, T, TSIZE, WORK, LWORK, */
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
/* > CGELQ computes a LQ factorization of an M-by-N matrix A. */
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
/* >          On exit, the elements on and below the diagonal of the array */
/* >          contain the M-by-min(M,N) lower trapezoidal matrix L */
/* >          (L is lower triangular if M <= N); */
/* >          the elements above the diagonal are used to store part of the */
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
/* > creating any LQ factorization algorithm they wish. The triangular */
/* > (trapezoidal) L has to be stored in the lower part of A. The lower part of A */
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
/* >                           CLASWLQ or CGELQT */
/* > */
/* >  Depending on the matrix dimensions M and N, and row and column */
/* >  block sizes MB and NB returned by ILAENV, CGELQ will use either */
/* >  CLASWLQ (if the matrix is short-and-wide) or CGELQT to compute */
/* >  the LQ factorization. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int cgelq_(integer *m, integer *n, doublecomplex *a, integer 
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
    extern /* Subroutine */ int cgelqt_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static logical lminws, lquery;
    static integer mintsz;
    extern /* Subroutine */ int claswlq_(integer *, integer *, integer *, 
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

#line 199 "cgelq.f"
    /* Parameter adjustments */
#line 199 "cgelq.f"
    a_dim1 = *lda;
#line 199 "cgelq.f"
    a_offset = 1 + a_dim1;
#line 199 "cgelq.f"
    a -= a_offset;
#line 199 "cgelq.f"
    --t;
#line 199 "cgelq.f"
    --work;
#line 199 "cgelq.f"

#line 199 "cgelq.f"
    /* Function Body */
#line 199 "cgelq.f"
    *info = 0;

#line 201 "cgelq.f"
    lquery = *tsize == -1 || *tsize == -2 || *lwork == -1 || *lwork == -2;

#line 204 "cgelq.f"
    mint = FALSE_;
#line 205 "cgelq.f"
    minw = FALSE_;
#line 206 "cgelq.f"
    if (*tsize == -2 || *lwork == -2) {
#line 207 "cgelq.f"
	if (*tsize != -1) {
#line 207 "cgelq.f"
	    mint = TRUE_;
#line 207 "cgelq.f"
	}
#line 208 "cgelq.f"
	if (*lwork != -1) {
#line 208 "cgelq.f"
	    minw = TRUE_;
#line 208 "cgelq.f"
	}
#line 209 "cgelq.f"
    }

/*     Determine the block size */

#line 213 "cgelq.f"
    if (min(*m,*n) > 0) {
#line 214 "cgelq.f"
	mb = ilaenv_(&c__1, "CGELQ ", " ", m, n, &c__1, &c_n1, (ftnlen)6, (
		ftnlen)1);
#line 215 "cgelq.f"
	nb = ilaenv_(&c__1, "CGELQ ", " ", m, n, &c__2, &c_n1, (ftnlen)6, (
		ftnlen)1);
#line 216 "cgelq.f"
    } else {
#line 217 "cgelq.f"
	mb = 1;
#line 218 "cgelq.f"
	nb = *n;
#line 219 "cgelq.f"
    }
#line 220 "cgelq.f"
    if (mb > min(*m,*n) || mb < 1) {
#line 220 "cgelq.f"
	mb = 1;
#line 220 "cgelq.f"
    }
#line 221 "cgelq.f"
    if (nb > *n || nb <= *m) {
#line 221 "cgelq.f"
	nb = *n;
#line 221 "cgelq.f"
    }
#line 222 "cgelq.f"
    mintsz = *m + 5;
#line 223 "cgelq.f"
    if (nb > *m && *n > *m) {
#line 224 "cgelq.f"
	if ((*n - *m) % (nb - *m) == 0) {
#line 225 "cgelq.f"
	    nblcks = (*n - *m) / (nb - *m);
#line 226 "cgelq.f"
	} else {
#line 227 "cgelq.f"
	    nblcks = (*n - *m) / (nb - *m) + 1;
#line 228 "cgelq.f"
	}
#line 229 "cgelq.f"
    } else {
#line 230 "cgelq.f"
	nblcks = 1;
#line 231 "cgelq.f"
    }

/*     Determine if the workspace size satisfies minimal size */

#line 235 "cgelq.f"
    lminws = FALSE_;
/* Computing MAX */
#line 236 "cgelq.f"
    i__1 = 1, i__2 = mb * *m * nblcks + 5;
#line 236 "cgelq.f"
    if ((*tsize < max(i__1,i__2) || *lwork < mb * *m) && *lwork >= *m && *
	    tsize >= mintsz && ! lquery) {
/* Computing MAX */
#line 239 "cgelq.f"
	i__1 = 1, i__2 = mb * *m * nblcks + 5;
#line 239 "cgelq.f"
	if (*tsize < max(i__1,i__2)) {
#line 240 "cgelq.f"
	    lminws = TRUE_;
#line 241 "cgelq.f"
	    mb = 1;
#line 242 "cgelq.f"
	    nb = *n;
#line 243 "cgelq.f"
	}
#line 244 "cgelq.f"
	if (*lwork < mb * *m) {
#line 245 "cgelq.f"
	    lminws = TRUE_;
#line 246 "cgelq.f"
	    mb = 1;
#line 247 "cgelq.f"
	}
#line 248 "cgelq.f"
    }

#line 250 "cgelq.f"
    if (*m < 0) {
#line 251 "cgelq.f"
	*info = -1;
#line 252 "cgelq.f"
    } else if (*n < 0) {
#line 253 "cgelq.f"
	*info = -2;
#line 254 "cgelq.f"
    } else if (*lda < max(1,*m)) {
#line 255 "cgelq.f"
	*info = -4;
#line 256 "cgelq.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 256 "cgelq.f"
	i__1 = 1, i__2 = mb * *m * nblcks + 5;
#line 256 "cgelq.f"
	if (*tsize < max(i__1,i__2) && ! lquery && ! lminws) {
#line 258 "cgelq.f"
	    *info = -6;
#line 259 "cgelq.f"
	} else /* if(complicated condition) */ {
/* Computing MAX */
#line 259 "cgelq.f"
	    i__1 = 1, i__2 = *m * mb;
#line 259 "cgelq.f"
	    if (*lwork < max(i__1,i__2) && ! lquery && ! lminws) {
#line 261 "cgelq.f"
		*info = -8;
#line 262 "cgelq.f"
	    }
#line 262 "cgelq.f"
	}
#line 262 "cgelq.f"
    }

#line 264 "cgelq.f"
    if (*info == 0) {
#line 265 "cgelq.f"
	if (mint) {
#line 266 "cgelq.f"
	    t[1].r = (doublereal) mintsz, t[1].i = 0.;
#line 267 "cgelq.f"
	} else {
#line 268 "cgelq.f"
	    i__1 = mb * *m * nblcks + 5;
#line 268 "cgelq.f"
	    t[1].r = (doublereal) i__1, t[1].i = 0.;
#line 269 "cgelq.f"
	}
#line 270 "cgelq.f"
	t[2].r = (doublereal) mb, t[2].i = 0.;
#line 271 "cgelq.f"
	t[3].r = (doublereal) nb, t[3].i = 0.;
#line 272 "cgelq.f"
	if (minw) {
#line 273 "cgelq.f"
	    i__1 = max(1,*n);
#line 273 "cgelq.f"
	    work[1].r = (doublereal) i__1, work[1].i = 0.;
#line 274 "cgelq.f"
	} else {
/* Computing MAX */
#line 275 "cgelq.f"
	    i__2 = 1, i__3 = mb * *m;
#line 275 "cgelq.f"
	    i__1 = max(i__2,i__3);
#line 275 "cgelq.f"
	    work[1].r = (doublereal) i__1, work[1].i = 0.;
#line 276 "cgelq.f"
	}
#line 277 "cgelq.f"
    }
#line 278 "cgelq.f"
    if (*info != 0) {
#line 279 "cgelq.f"
	i__1 = -(*info);
#line 279 "cgelq.f"
	xerbla_("CGELQ", &i__1, (ftnlen)5);
#line 280 "cgelq.f"
	return 0;
#line 281 "cgelq.f"
    } else if (lquery) {
#line 282 "cgelq.f"
	return 0;
#line 283 "cgelq.f"
    }

/*     Quick return if possible */

#line 287 "cgelq.f"
    if (min(*m,*n) == 0) {
#line 288 "cgelq.f"
	return 0;
#line 289 "cgelq.f"
    }

/*     The LQ Decomposition */

#line 293 "cgelq.f"
    if (*n <= *m || nb <= *m || nb >= *n) {
#line 294 "cgelq.f"
	cgelqt_(m, n, &mb, &a[a_offset], lda, &t[6], &mb, &work[1], info);
#line 295 "cgelq.f"
    } else {
#line 296 "cgelq.f"
	claswlq_(m, n, &mb, &nb, &a[a_offset], lda, &t[6], &mb, &work[1], 
		lwork, info);
#line 298 "cgelq.f"
    }

/* Computing MAX */
#line 300 "cgelq.f"
    i__2 = 1, i__3 = mb * *m;
#line 300 "cgelq.f"
    i__1 = max(i__2,i__3);
#line 300 "cgelq.f"
    work[1].r = (doublereal) i__1, work[1].i = 0.;

#line 302 "cgelq.f"
    return 0;

/*     End of CGELQ */

} /* cgelq_ */

