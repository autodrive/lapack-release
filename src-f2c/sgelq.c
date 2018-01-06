#line 1 "sgelq.f"
/* sgelq.f -- translated by f2c (version 20100827).
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

#line 1 "sgelq.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;


/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGELQ( M, N, A, LDA, T, TSIZE, WORK, LWORK, */
/*                         INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER           INFO, LDA, M, N, TSIZE, LWORK */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL              A( LDA, * ), T( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > SGELQ computes a LQ factorization of an M-by-N matrix A. */
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
/* >          A is REAL array, dimension (LDA,N) */
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
/* >          T is REAL array, dimension (MAX(5,TSIZE)) */
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
/* >          (workspace) REAL array, dimension (MAX(1,LWORK)) */
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
/* >                           SLASWLQ or SGELQT */
/* > */
/* >  Depending on the matrix dimensions M and N, and row and column */
/* >  block sizes MB and NB returned by ILAENV, SGELQ will use either */
/* >  SLASWLQ (if the matrix is short-and-wide) or SGELQT to compute */
/* >  the LQ factorization. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int sgelq_(integer *m, integer *n, doublereal *a, integer *
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
    extern /* Subroutine */ int sgelqt_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static logical lminws, lquery;
    static integer mintsz;
    extern /* Subroutine */ int slaswlq_(integer *, integer *, integer *, 
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
/*     .. Executable statements .. */

/*     Test the input arguments */

#line 199 "sgelq.f"
    /* Parameter adjustments */
#line 199 "sgelq.f"
    a_dim1 = *lda;
#line 199 "sgelq.f"
    a_offset = 1 + a_dim1;
#line 199 "sgelq.f"
    a -= a_offset;
#line 199 "sgelq.f"
    --t;
#line 199 "sgelq.f"
    --work;
#line 199 "sgelq.f"

#line 199 "sgelq.f"
    /* Function Body */
#line 199 "sgelq.f"
    *info = 0;

#line 201 "sgelq.f"
    lquery = *tsize == -1 || *tsize == -2 || *lwork == -1 || *lwork == -2;

#line 204 "sgelq.f"
    mint = FALSE_;
#line 205 "sgelq.f"
    minw = FALSE_;
#line 206 "sgelq.f"
    if (*tsize == -2 || *lwork == -2) {
#line 207 "sgelq.f"
	if (*tsize != -1) {
#line 207 "sgelq.f"
	    mint = TRUE_;
#line 207 "sgelq.f"
	}
#line 208 "sgelq.f"
	if (*lwork != -1) {
#line 208 "sgelq.f"
	    minw = TRUE_;
#line 208 "sgelq.f"
	}
#line 209 "sgelq.f"
    }

/*     Determine the block size */

#line 213 "sgelq.f"
    if (min(*m,*n) > 0) {
#line 214 "sgelq.f"
	mb = ilaenv_(&c__1, "SGELQ ", " ", m, n, &c__1, &c_n1, (ftnlen)6, (
		ftnlen)1);
#line 215 "sgelq.f"
	nb = ilaenv_(&c__1, "SGELQ ", " ", m, n, &c__2, &c_n1, (ftnlen)6, (
		ftnlen)1);
#line 216 "sgelq.f"
    } else {
#line 217 "sgelq.f"
	mb = 1;
#line 218 "sgelq.f"
	nb = *n;
#line 219 "sgelq.f"
    }
#line 220 "sgelq.f"
    if (mb > min(*m,*n) || mb < 1) {
#line 220 "sgelq.f"
	mb = 1;
#line 220 "sgelq.f"
    }
#line 221 "sgelq.f"
    if (nb > *n || nb <= *m) {
#line 221 "sgelq.f"
	nb = *n;
#line 221 "sgelq.f"
    }
#line 222 "sgelq.f"
    mintsz = *m + 5;
#line 223 "sgelq.f"
    if (nb > *m && *n > *m) {
#line 224 "sgelq.f"
	if ((*n - *m) % (nb - *m) == 0) {
#line 225 "sgelq.f"
	    nblcks = (*n - *m) / (nb - *m);
#line 226 "sgelq.f"
	} else {
#line 227 "sgelq.f"
	    nblcks = (*n - *m) / (nb - *m) + 1;
#line 228 "sgelq.f"
	}
#line 229 "sgelq.f"
    } else {
#line 230 "sgelq.f"
	nblcks = 1;
#line 231 "sgelq.f"
    }

/*     Determine if the workspace size satisfies minimal size */

#line 235 "sgelq.f"
    lminws = FALSE_;
/* Computing MAX */
#line 236 "sgelq.f"
    i__1 = 1, i__2 = mb * *m * nblcks + 5;
#line 236 "sgelq.f"
    if ((*tsize < max(i__1,i__2) || *lwork < mb * *m) && *lwork >= *m && *
	    tsize >= mintsz && ! lquery) {
/* Computing MAX */
#line 239 "sgelq.f"
	i__1 = 1, i__2 = mb * *m * nblcks + 5;
#line 239 "sgelq.f"
	if (*tsize < max(i__1,i__2)) {
#line 240 "sgelq.f"
	    lminws = TRUE_;
#line 241 "sgelq.f"
	    mb = 1;
#line 242 "sgelq.f"
	    nb = *n;
#line 243 "sgelq.f"
	}
#line 244 "sgelq.f"
	if (*lwork < mb * *m) {
#line 245 "sgelq.f"
	    lminws = TRUE_;
#line 246 "sgelq.f"
	    mb = 1;
#line 247 "sgelq.f"
	}
#line 248 "sgelq.f"
    }

#line 250 "sgelq.f"
    if (*m < 0) {
#line 251 "sgelq.f"
	*info = -1;
#line 252 "sgelq.f"
    } else if (*n < 0) {
#line 253 "sgelq.f"
	*info = -2;
#line 254 "sgelq.f"
    } else if (*lda < max(1,*m)) {
#line 255 "sgelq.f"
	*info = -4;
#line 256 "sgelq.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 256 "sgelq.f"
	i__1 = 1, i__2 = mb * *m * nblcks + 5;
#line 256 "sgelq.f"
	if (*tsize < max(i__1,i__2) && ! lquery && ! lminws) {
#line 258 "sgelq.f"
	    *info = -6;
#line 259 "sgelq.f"
	} else /* if(complicated condition) */ {
/* Computing MAX */
#line 259 "sgelq.f"
	    i__1 = 1, i__2 = *m * mb;
#line 259 "sgelq.f"
	    if (*lwork < max(i__1,i__2) && ! lquery && ! lminws) {
#line 261 "sgelq.f"
		*info = -8;
#line 262 "sgelq.f"
	    }
#line 262 "sgelq.f"
	}
#line 262 "sgelq.f"
    }

#line 264 "sgelq.f"
    if (*info == 0) {
#line 265 "sgelq.f"
	if (mint) {
#line 266 "sgelq.f"
	    t[1] = (doublereal) mintsz;
#line 267 "sgelq.f"
	} else {
#line 268 "sgelq.f"
	    t[1] = (doublereal) (mb * *m * nblcks + 5);
#line 269 "sgelq.f"
	}
#line 270 "sgelq.f"
	t[2] = (doublereal) mb;
#line 271 "sgelq.f"
	t[3] = (doublereal) nb;
#line 272 "sgelq.f"
	if (minw) {
#line 273 "sgelq.f"
	    work[1] = (doublereal) max(1,*n);
#line 274 "sgelq.f"
	} else {
/* Computing MAX */
#line 275 "sgelq.f"
	    i__1 = 1, i__2 = mb * *m;
#line 275 "sgelq.f"
	    work[1] = (doublereal) max(i__1,i__2);
#line 276 "sgelq.f"
	}
#line 277 "sgelq.f"
    }
#line 278 "sgelq.f"
    if (*info != 0) {
#line 279 "sgelq.f"
	i__1 = -(*info);
#line 279 "sgelq.f"
	xerbla_("SGELQ", &i__1, (ftnlen)5);
#line 280 "sgelq.f"
	return 0;
#line 281 "sgelq.f"
    } else if (lquery) {
#line 282 "sgelq.f"
	return 0;
#line 283 "sgelq.f"
    }

/*     Quick return if possible */

#line 287 "sgelq.f"
    if (min(*m,*n) == 0) {
#line 288 "sgelq.f"
	return 0;
#line 289 "sgelq.f"
    }

/*     The LQ Decomposition */

#line 293 "sgelq.f"
    if (*n <= *m || nb <= *m || nb >= *n) {
#line 294 "sgelq.f"
	sgelqt_(m, n, &mb, &a[a_offset], lda, &t[6], &mb, &work[1], info);
#line 295 "sgelq.f"
    } else {
#line 296 "sgelq.f"
	slaswlq_(m, n, &mb, &nb, &a[a_offset], lda, &t[6], &mb, &work[1], 
		lwork, info);
#line 298 "sgelq.f"
    }

/* Computing MAX */
#line 300 "sgelq.f"
    i__1 = 1, i__2 = mb * *m;
#line 300 "sgelq.f"
    work[1] = (doublereal) max(i__1,i__2);
#line 301 "sgelq.f"
    return 0;

/*     End of SGELQ */

} /* sgelq_ */

