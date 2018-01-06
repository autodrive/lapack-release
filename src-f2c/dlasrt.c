#line 1 "dlasrt.f"
/* dlasrt.f -- translated by f2c (version 20100827).
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

#line 1 "dlasrt.f"
/* > \brief \b DLASRT sorts numbers in increasing or decreasing order. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLASRT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasrt.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasrt.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasrt.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLASRT( ID, N, D, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          ID */
/*       INTEGER            INFO, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   D( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > Sort the numbers in D in increasing order (if ID = 'I') or */
/* > in decreasing order (if ID = 'D' ). */
/* > */
/* > Use Quick Sort, reverting to Insertion sort on arrays of */
/* > size <= 20. Dimension of STACK limits N to about 2**32. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] ID */
/* > \verbatim */
/* >          ID is CHARACTER*1 */
/* >          = 'I': sort D in increasing order; */
/* >          = 'D': sort D in decreasing order. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The length of the array D. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          On entry, the array to be sorted. */
/* >          On exit, D has been sorted into increasing order */
/* >          (D(1) <= ... <= D(N) ) or into decreasing order */
/* >          (D(1) >= ... >= D(N) ), depending on ID. */
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

/* > \date September 2012 */

/* > \ingroup auxOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int dlasrt_(char *id, integer *n, doublereal *d__, integer *
	info, ftnlen id_len)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j;
    static doublereal d1, d2, d3;
    static integer dir;
    static doublereal tmp;
    static integer endd;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer stack[64]	/* was [2][32] */;
    static doublereal dmnmx;
    static integer start;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer stkpnt;


/*  -- LAPACK computational routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input paramters. */

#line 128 "dlasrt.f"
    /* Parameter adjustments */
#line 128 "dlasrt.f"
    --d__;
#line 128 "dlasrt.f"

#line 128 "dlasrt.f"
    /* Function Body */
#line 128 "dlasrt.f"
    *info = 0;
#line 129 "dlasrt.f"
    dir = -1;
#line 130 "dlasrt.f"
    if (lsame_(id, "D", (ftnlen)1, (ftnlen)1)) {
#line 131 "dlasrt.f"
	dir = 0;
#line 132 "dlasrt.f"
    } else if (lsame_(id, "I", (ftnlen)1, (ftnlen)1)) {
#line 133 "dlasrt.f"
	dir = 1;
#line 134 "dlasrt.f"
    }
#line 135 "dlasrt.f"
    if (dir == -1) {
#line 136 "dlasrt.f"
	*info = -1;
#line 137 "dlasrt.f"
    } else if (*n < 0) {
#line 138 "dlasrt.f"
	*info = -2;
#line 139 "dlasrt.f"
    }
#line 140 "dlasrt.f"
    if (*info != 0) {
#line 141 "dlasrt.f"
	i__1 = -(*info);
#line 141 "dlasrt.f"
	xerbla_("DLASRT", &i__1, (ftnlen)6);
#line 142 "dlasrt.f"
	return 0;
#line 143 "dlasrt.f"
    }

/*     Quick return if possible */

#line 147 "dlasrt.f"
    if (*n <= 1) {
#line 147 "dlasrt.f"
	return 0;
#line 147 "dlasrt.f"
    }

#line 150 "dlasrt.f"
    stkpnt = 1;
#line 151 "dlasrt.f"
    stack[0] = 1;
#line 152 "dlasrt.f"
    stack[1] = *n;
#line 153 "dlasrt.f"
L10:
#line 154 "dlasrt.f"
    start = stack[(stkpnt << 1) - 2];
#line 155 "dlasrt.f"
    endd = stack[(stkpnt << 1) - 1];
#line 156 "dlasrt.f"
    --stkpnt;
#line 157 "dlasrt.f"
    if (endd - start <= 20 && endd - start > 0) {

/*        Do Insertion sort on D( START:ENDD ) */

#line 161 "dlasrt.f"
	if (dir == 0) {

/*           Sort into decreasing order */

#line 165 "dlasrt.f"
	    i__1 = endd;
#line 165 "dlasrt.f"
	    for (i__ = start + 1; i__ <= i__1; ++i__) {
#line 166 "dlasrt.f"
		i__2 = start + 1;
#line 166 "dlasrt.f"
		for (j = i__; j >= i__2; --j) {
#line 167 "dlasrt.f"
		    if (d__[j] > d__[j - 1]) {
#line 168 "dlasrt.f"
			dmnmx = d__[j];
#line 169 "dlasrt.f"
			d__[j] = d__[j - 1];
#line 170 "dlasrt.f"
			d__[j - 1] = dmnmx;
#line 171 "dlasrt.f"
		    } else {
#line 172 "dlasrt.f"
			goto L30;
#line 173 "dlasrt.f"
		    }
#line 174 "dlasrt.f"
/* L20: */
#line 174 "dlasrt.f"
		}
#line 175 "dlasrt.f"
L30:
#line 175 "dlasrt.f"
		;
#line 175 "dlasrt.f"
	    }

#line 177 "dlasrt.f"
	} else {

/*           Sort into increasing order */

#line 181 "dlasrt.f"
	    i__1 = endd;
#line 181 "dlasrt.f"
	    for (i__ = start + 1; i__ <= i__1; ++i__) {
#line 182 "dlasrt.f"
		i__2 = start + 1;
#line 182 "dlasrt.f"
		for (j = i__; j >= i__2; --j) {
#line 183 "dlasrt.f"
		    if (d__[j] < d__[j - 1]) {
#line 184 "dlasrt.f"
			dmnmx = d__[j];
#line 185 "dlasrt.f"
			d__[j] = d__[j - 1];
#line 186 "dlasrt.f"
			d__[j - 1] = dmnmx;
#line 187 "dlasrt.f"
		    } else {
#line 188 "dlasrt.f"
			goto L50;
#line 189 "dlasrt.f"
		    }
#line 190 "dlasrt.f"
/* L40: */
#line 190 "dlasrt.f"
		}
#line 191 "dlasrt.f"
L50:
#line 191 "dlasrt.f"
		;
#line 191 "dlasrt.f"
	    }

#line 193 "dlasrt.f"
	}

#line 195 "dlasrt.f"
    } else if (endd - start > 20) {

/*        Partition D( START:ENDD ) and stack parts, largest one first */

/*        Choose partition entry as median of 3 */

#line 201 "dlasrt.f"
	d1 = d__[start];
#line 202 "dlasrt.f"
	d2 = d__[endd];
#line 203 "dlasrt.f"
	i__ = (start + endd) / 2;
#line 204 "dlasrt.f"
	d3 = d__[i__];
#line 205 "dlasrt.f"
	if (d1 < d2) {
#line 206 "dlasrt.f"
	    if (d3 < d1) {
#line 207 "dlasrt.f"
		dmnmx = d1;
#line 208 "dlasrt.f"
	    } else if (d3 < d2) {
#line 209 "dlasrt.f"
		dmnmx = d3;
#line 210 "dlasrt.f"
	    } else {
#line 211 "dlasrt.f"
		dmnmx = d2;
#line 212 "dlasrt.f"
	    }
#line 213 "dlasrt.f"
	} else {
#line 214 "dlasrt.f"
	    if (d3 < d2) {
#line 215 "dlasrt.f"
		dmnmx = d2;
#line 216 "dlasrt.f"
	    } else if (d3 < d1) {
#line 217 "dlasrt.f"
		dmnmx = d3;
#line 218 "dlasrt.f"
	    } else {
#line 219 "dlasrt.f"
		dmnmx = d1;
#line 220 "dlasrt.f"
	    }
#line 221 "dlasrt.f"
	}

#line 223 "dlasrt.f"
	if (dir == 0) {

/*           Sort into decreasing order */

#line 227 "dlasrt.f"
	    i__ = start - 1;
#line 228 "dlasrt.f"
	    j = endd + 1;
#line 229 "dlasrt.f"
L60:
#line 230 "dlasrt.f"
L70:
#line 231 "dlasrt.f"
	    --j;
#line 232 "dlasrt.f"
	    if (d__[j] < dmnmx) {
#line 232 "dlasrt.f"
		goto L70;
#line 232 "dlasrt.f"
	    }
#line 234 "dlasrt.f"
L80:
#line 235 "dlasrt.f"
	    ++i__;
#line 236 "dlasrt.f"
	    if (d__[i__] > dmnmx) {
#line 236 "dlasrt.f"
		goto L80;
#line 236 "dlasrt.f"
	    }
#line 238 "dlasrt.f"
	    if (i__ < j) {
#line 239 "dlasrt.f"
		tmp = d__[i__];
#line 240 "dlasrt.f"
		d__[i__] = d__[j];
#line 241 "dlasrt.f"
		d__[j] = tmp;
#line 242 "dlasrt.f"
		goto L60;
#line 243 "dlasrt.f"
	    }
#line 244 "dlasrt.f"
	    if (j - start > endd - j - 1) {
#line 245 "dlasrt.f"
		++stkpnt;
#line 246 "dlasrt.f"
		stack[(stkpnt << 1) - 2] = start;
#line 247 "dlasrt.f"
		stack[(stkpnt << 1) - 1] = j;
#line 248 "dlasrt.f"
		++stkpnt;
#line 249 "dlasrt.f"
		stack[(stkpnt << 1) - 2] = j + 1;
#line 250 "dlasrt.f"
		stack[(stkpnt << 1) - 1] = endd;
#line 251 "dlasrt.f"
	    } else {
#line 252 "dlasrt.f"
		++stkpnt;
#line 253 "dlasrt.f"
		stack[(stkpnt << 1) - 2] = j + 1;
#line 254 "dlasrt.f"
		stack[(stkpnt << 1) - 1] = endd;
#line 255 "dlasrt.f"
		++stkpnt;
#line 256 "dlasrt.f"
		stack[(stkpnt << 1) - 2] = start;
#line 257 "dlasrt.f"
		stack[(stkpnt << 1) - 1] = j;
#line 258 "dlasrt.f"
	    }
#line 259 "dlasrt.f"
	} else {

/*           Sort into increasing order */

#line 263 "dlasrt.f"
	    i__ = start - 1;
#line 264 "dlasrt.f"
	    j = endd + 1;
#line 265 "dlasrt.f"
L90:
#line 266 "dlasrt.f"
L100:
#line 267 "dlasrt.f"
	    --j;
#line 268 "dlasrt.f"
	    if (d__[j] > dmnmx) {
#line 268 "dlasrt.f"
		goto L100;
#line 268 "dlasrt.f"
	    }
#line 270 "dlasrt.f"
L110:
#line 271 "dlasrt.f"
	    ++i__;
#line 272 "dlasrt.f"
	    if (d__[i__] < dmnmx) {
#line 272 "dlasrt.f"
		goto L110;
#line 272 "dlasrt.f"
	    }
#line 274 "dlasrt.f"
	    if (i__ < j) {
#line 275 "dlasrt.f"
		tmp = d__[i__];
#line 276 "dlasrt.f"
		d__[i__] = d__[j];
#line 277 "dlasrt.f"
		d__[j] = tmp;
#line 278 "dlasrt.f"
		goto L90;
#line 279 "dlasrt.f"
	    }
#line 280 "dlasrt.f"
	    if (j - start > endd - j - 1) {
#line 281 "dlasrt.f"
		++stkpnt;
#line 282 "dlasrt.f"
		stack[(stkpnt << 1) - 2] = start;
#line 283 "dlasrt.f"
		stack[(stkpnt << 1) - 1] = j;
#line 284 "dlasrt.f"
		++stkpnt;
#line 285 "dlasrt.f"
		stack[(stkpnt << 1) - 2] = j + 1;
#line 286 "dlasrt.f"
		stack[(stkpnt << 1) - 1] = endd;
#line 287 "dlasrt.f"
	    } else {
#line 288 "dlasrt.f"
		++stkpnt;
#line 289 "dlasrt.f"
		stack[(stkpnt << 1) - 2] = j + 1;
#line 290 "dlasrt.f"
		stack[(stkpnt << 1) - 1] = endd;
#line 291 "dlasrt.f"
		++stkpnt;
#line 292 "dlasrt.f"
		stack[(stkpnt << 1) - 2] = start;
#line 293 "dlasrt.f"
		stack[(stkpnt << 1) - 1] = j;
#line 294 "dlasrt.f"
	    }
#line 295 "dlasrt.f"
	}
#line 296 "dlasrt.f"
    }
#line 297 "dlasrt.f"
    if (stkpnt > 0) {
#line 297 "dlasrt.f"
	goto L10;
#line 297 "dlasrt.f"
    }
#line 299 "dlasrt.f"
    return 0;

/*     End of DLASRT */

} /* dlasrt_ */

