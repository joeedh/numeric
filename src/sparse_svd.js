let abs   = Math.abs,
    fabs  = Math.abs,
    sqrt  = Math.sqrt,
    floor = Math.floor,
    ceil  = Math.ceil,
    atan  = Math.atan,
    pow   = Math.pow,
    log   = Math.log;
let eps, eps1, reps, eps34;
let LanStore, OPBTemp;

const MAXLL = 2
const LMTNW = 100000000 /* max. size of working area allowed  */

let a = 0;
const STORQ = a++, RETRQ = a++, STORP = a++, RETRP = a++;

export var ierr = 0;

import {
  SMat, DMat, SVDRec, cachering, svd_datx, svd_daxpy, svd_dcopy,
  svd_ddot, svd_dscal, svd_dsort2, svd_dswap, svd_fsign, svd_idamax, svd_opa,
  svd_opb, svd_pythag, svdWriteDenseArray, SvdError, ArrayPool, dtempArray,
  itempArray, SAFE_FREE, svd_doubleArray, svd_longArray, tempArray
} from './sparse_svd_util.js';

export * from './sparse_svd_util.js';


function array2(a, b) {
  let ret = tempArray(2);
  ret[0] = a;
  ret[1] = b;
  return ret;
}

function array3(a, b, c) {
  let ret = tempArray(3);
  ret[0] = a;
  ret[1] = b;
  ret[2] = c;
  return ret;
}


/***********************************************************************
 *                                                                     *
 *        random()                               *
 *                        (double precision)                           *
 ***********************************************************************/
/***********************************************************************

 Description
 -----------

 This is a translation of a Fortran-77 uniform random number
 generator.  The code is based  on  theory and suggestions  given in
 D. E. Knuth (1969),  vol  2.  The argument to the function should
 be initialized to an arbitrary integer prior to the first call to
 random.  The calling program should  not  alter  the value of the
 argument between subsequent calls to random.  Random returns values
 within the interval (0,1).


 Arguments
 ---------

 (input)
 iy     an integer seed whose value must not be altered by the caller
 between subsequent calls

 (output)
 random  a double precision random number between (0,1)

 ***********************************************************************/
export const svd_random2 = (function () {
  let m2 = 0, ia = 0, ic = 0, mic = 0, halfm = 0, s = 0;

  return function svd_random2(iy) {
    /* If first entry, compute (max int) / 2 */
    if (!m2) {
      m2 = 1<<(8*4 - 2);
      halfm = m2;

      /* compute multiplier and increment for linear congruential
       * method */
      ia = 8* ~~(halfm*atan(1.0)/8.0) + 5;
      ic = 2* ~~(halfm*(0.5 - sqrt(3.0)/6.0)) + 1;
      mic = (m2 - ic) + m2;

      /* s is the scale factor for converting to floating point */
      s = 0.5/halfm;
    }

    /* compute next random number */
    iy.i = iy.i*ia;

    /* for computers which do not allow integer overflow on addition */
    if (iy.i > mic) iy.i = (iy.i - m2) - m2;

    iy.i = iy.i + ic;

    /* for computers whose word length for addition is greater than
     * for multiplication */
    if (iy.i/2 > m2) iy.i = (iy.i - m2) - m2;

    /* for computers whose integer overflow affects the sign bit */
    if (iy.i < 0) iy.i = (iy.i + m2) + m2;

    return (iy.i*s);
  }
})();

//SVDRec svdLAS2A(SMat A, long dimensions) {
let _etmp = [0, 0];

export function svdLAS2A(A, dimensions) {
  let end = _etmp;

  end[0] = -1.0e-30;
  end[1] = 1.0e-30;

  let kappa = 1e-6;

  return svdLAS2(A, dimensions, 0, end, kappa);
}

export var SVDVerbosity = 0;

function printf() {
  console.log(...arguments);
}

/**************************************************************
 * Function copies a vector x to a vector y              *
 * Based on Fortran-77 routine from Linpack by J. Dongarra    *
 **************************************************************/
//void svd_dcopy(long n, double *dx, long incx, double *dy, long incy) {
function svd_dcopy(n, dx, incx, dy, incy, startx = 0, starty = 0) {
  long
  i;
  if (!incx || !incy || n <= 0) {
    return;
  }

  if (incx === 1 && incy === 1) {
    let j = starty;
    for (let i = startx; i < startx + n; i++, j++) {
      dx[i] = dy[j];
    }

    return;
  }

  let j = starty;
  for (let i = startx; i < startx + n; i += incx, j += incy) {
    dx[i] = dy[j];
  }
}


function memset(arr, val, len) {
  if (arr instanceof Uint8Array || arr instanceof ClampedUint8Array) {
    for (let i = 0; i < len; i++) {
      arr[i] = val;
    }
  } else if (arr.buffer && !(arr instanceof Array)) {
    arr = new Uint8Array(arr.buffer);
    return memset(arr, val, len);
  } else {
    throw new Error("memset not called with typed array");
  }
}

//SVDRec svdLAS2(SMat A, long dimensions, long iterations, double end[2],
//               double kappa) {
export function svdLAS2(A, dimensions, iterations, end, kappa) {
  let tranpose = false;
  let n, i, steps, nsig, neig, m;
  let wptr = tempArray(10), ritz, bnd;
  let R;

  m = Math.min(A.rows, A.cols);

  if (dimensions <= 0 || dimensions > m)
    dimensions = m;
  if (iterations <= 0 || iterations > m)
    iterations = m;
  if (iterations < dimensions) iterations = dimensions;

  /* Check parameters */
  if (check_parameters(A, dimensions, iterations, end[0], end[1], true))
    return undefined;

  let transpose = false;

  /* If A is wide, the SVD is computed on its transpose for speed. */
  if (A.cols >= A.rows*1.2) {
    if (SVDVerbosity > 0) printf("TRANSPOSING THE MATRIX FOR SPEED\n");
    transpose = true;
    A = svdTransposeS(A);
  }

  n = A.cols;
  /* Compute machine precision */
  let {ibeta, it, irnd, machep, negep, eps0} = machar();

  eps = eps0;
  eps1 = eps*sqrt(n);
  reps = sqrt(eps);
  eps34 = reps*sqrt(reps);

  /* Allocate temporary space. */
  wptr[0] = svd_doubleArray(n, true, "las2: wptr[0]");
  wptr[1] = svd_doubleArray(n, false, "las2: wptr[1]");
  wptr[2] = svd_doubleArray(n, false, "las2: wptr[2]");
  wptr[3] = svd_doubleArray(n, false, "las2: wptr[3]");
  wptr[4] = svd_doubleArray(n, false, "las2: wptr[4]");
  wptr[5] = svd_doubleArray(n, false, "las2: wptr[5]");
  wptr[6] = svd_doubleArray(iterations, false, "las2: wptr[6]");

  wptr[7] = svd_doubleArray(iterations, false, "las2: wptr[7]");

  wptr[8] = svd_doubleArray(iterations, false, "las2: wptr[8]");

  wptr[9] = svd_doubleArray(iterations + 1, false, "las2: wptr[9]");

  /* Calloc may be unnecessary: */
  ritz = svd_doubleArray(iterations + 1, true, "las2: ritz");

  /* Calloc may be unnecessary: */
  bnd = svd_doubleArray(iterations + 1, true, "las2: bnd");

  memset(bnd, 127, (iterations + 1)*8);

  LanStore = tempArray(iterations + MAXLL);
  OPBTemp = svd_doubleArray(A.rows, false, "las2: OPBTemp");

  ierr = 0;

  /* Actually run the lanczos thing: */
  let ret = lanso(A, iterations, dimensions, end[0], end[1], ritz, bnd, wptr, n);
  steps = ret.steps;
  neig = ret.neig;

  /* Print some stuff. */
  if (SVDVerbosity > 0) {
    printf(`
  NUMBER OF LANCZOS STEPS   = ${steps + 1}
  RITZ VALUES STABILIZED    = ${neig.toFixed(5)}\n`.trim()
    );
  }
  if (SVDVerbosity > 2) {
    printf("\nCOMPUTED RITZ VALUES  (ERROR BNDS)\n");
    for (let i = 0; i <= steps; i++) {
      printf(`${i + 1} ${ritz[i]} ${bnd[i]}\n`);
    }
  }

  function cleanup() {
    for (i = 0; i <= 9; i++) {
      SAFE_FREE(wptr[i]);
    }
    SAFE_FREE(ritz);
    SAFE_FREE(bnd);
    if (LanStore) {
      for (i = 0; i < iterations + MAXLL; i++) {
        SAFE_FREE(LanStore[i]);
      }
      SAFE_FREE(LanStore);
    }
    SAFE_FREE(OPBTemp);
  }

  SAFE_FREE(wptr[0]);
  SAFE_FREE(wptr[1]);
  SAFE_FREE(wptr[2]);
  SAFE_FREE(wptr[3]);
  SAFE_FREE(wptr[4]);
  SAFE_FREE(wptr[7]);
  SAFE_FREE(wptr[8]);

  /* Compute eigenvectors */
  kappa = Math.max(fabs(kappa), eps34);

  R = svdNewSVDRec();
  if (!R) {
    cleanup();
    throw new SvdError("svdLAS2: allocation of R failed");
  }

  R.d = /*Math.min(nsig, dimensions)*/dimensions;
  R.Ut = svdNewDMat(R.d, A.rows);
  R.S = svd_doubleArray(R.d, true, "las2: R.s");
  R.Vt = svdNewDMat(R.d, A.cols);

  if (!R.Ut || !R.S || !R.Vt) {
    cleanup();
    throw new SvdError("svdLAS2: allocation of R failed");
  }

  nsig = ritvec(n, A, R, kappa, ritz, bnd, wptr[6], wptr[9], wptr[5], steps,
    neig);

  if (SVDVerbosity > 1) {
    printf("\nSINGULAR VALUES: ");
    svdWriteDenseArray(R.S, R.d, "-", false);

    if (SVDVerbosity > 2) {
      printf("\nLEFT SINGULAR VECTORS (transpose of U): ");
      svdWriteDenseMatrix(R.Ut, "-", SVD_F_DT);

      printf("\nRIGHT SINGULAR VECTORS (transpose of V): ");
      svdWriteDenseMatrix(R.Vt, "-", SVD_F_DT);
    }
    printf("SINGULAR VALUES FOUND     = %6d; nsig=%ld\n", R.d, nsig);
  } else if (SVDVerbosity > 0)
    printf("SINGULAR VALUES FOUND     = %6d\n", R.d);


  cleanup();

  /* This swaps and transposes the singular matrices if A was transposed. */
  if (R && transpose) {
    let T = new DMat();
    svdFreeSMat(A);
    T = R.Ut;
    R.Ut = R.Vt;
    R.Vt = T;
  }

  return R;
}


/***********************************************************************
 *                                                                     *
 *                        ritvec()                                     *
 *      Function computes the singular vectors of matrix A         *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

 Description
 -----------

 This function is invoked by landr() only if eigenvectors of the A'A
 eigenproblem are desired.  When called, ritvec() computes the
 singular vectors of A and writes the result to an unformatted file.


 Parameters
 ----------

 (input)
 nrow       number of rows of A
 steps      number of Lanczos iterations performed
 fp_out2    pointer to unformatted output file
 n        dimension of matrix A
 kappa      relative accuracy of ritz values acceptable as
 eigenvalues of A'A
 ritz       array of ritz values
 bnd        array of error bounds
 alf        array of diagonal elements of the tridiagonal matrix T
 bet        array of off-diagonal elements of T
 w1, w2     work space

 (output)
 xv1        array of eigenvectors of A'A (right singular vectors of A)
 ierr        error code
 0 for normal return from imtql2()
 k if convergence did not occur for k-th eigenvalue in
 imtql2()
 nsig       number of accepted ritz values based on kappa

 (local)
 s        work array which is initialized to the identity matrix
 of order (j + 1) upon calling imtql2().  After the call,
 s contains the orthonormal eigenvectors of the symmetric
 tridiagonal matrix T

 Functions used
 --------------

 BLAS    svd_dscal, svd_dcopy, svd_daxpy
 USER    store
 imtql2

 ***********************************************************************/

//void rotateArray(double *a, int size, int x) {
function rotateArray(a, size, x) {
  let i, j, n, start;
  let t1, t2;
  if (x === 0) return;
  j = start = 0;
  t1 = a[0];
  for (i = 0; i < size; i++) {
    n = (j >= x) ? j - x : j + size - x;
    t2 = a[n];
    a[n] = t1;
    t1 = t2;
    j = n;
    if (j === start) {
      start = ++j;
      t1 = a[j];
    }
  }
}

function ritvec(n, A, R, kappa, ritz, bnd, alf, bet, w2, steps, neig) {
//long ritvec(long n, SMat A, SVDRec R, double kappa, double *ritz, double *bnd,
//            double *alf, double *bet, double *w2, long steps, long neig) {
  let js, jsq, i, k, /*size,*/ id2, tmp, nsig, x;
  let s, xv2, tmp0, tmp1, xnorm, w1 = R.Vt.value[0];

  js = steps + 1;
  jsq = js*js;
  /*size = sizeof(double) * n;*/

  s = svd_doubleArray(jsq, true, "ritvec: s");
  xv2 = svd_doubleArray(n, false, "ritvec: xv2");

  /* initialize s to an identity matrix */
  for (i = 0; i < jsq; i += (js + 1)) {
    s[i] = 1.0;
  }

  svd_dcopy(js, alf, 1, w1, -1);
  svd_dcopy(steps, bet, 1, w2, -1, 1, 1);

  /* on return from imtql2(), w1 contains eigenvalues in ascending
   * order and s contains the corresponding eigenvectors */
  imtql2(js, js, w1, w2, s);

  /*fwrite((char *)&n, sizeof(n), 1, fp_out2);
    fwrite((char *)&js, sizeof(js), 1, fp_out2);
    fwrite((char *)&kappa, sizeof(kappa), 1, fp_out2);*/
  /*id = 0;*/
  nsig = 0;

  if (ierr) {
    R.d = 0;
  } else {
    x = 0;
    id2 = jsq - js;
    for (k = 0; k < js; k++) {
      tmp = id2;
      if (bnd[k] <= kappa*fabs(ritz[k]) && k > js - neig - 1) {
        if (--x < 0) x = R.d - 1;
        w1 = R.Vt.value[x];
        for (i = 0; i < n; i++) {
          w1[i] = 0.0;
        }
        for (i = 0; i < js; i++) {
          store(n, RETRQ, i, w2);
          svd_daxpy(n, s[tmp], w2, 1, w1, 1);
          tmp -= js;
        }
        /*fwrite((char *)w1, size, 1, fp_out2);*/

        /* store the w1 vector row-wise in array xv1;
         * size of xv1 is (steps+1) * (nrow+ncol) elements
         * and each vector, even though only ncol long,
         * will have (nrow+ncol) elements in xv1.
         * It is as if xv1 is a 2-d array (steps+1) by
         * (nrow+ncol) and each vector occupies a row  */

        /* j is the index in the R arrays, which are sorted by high to low
           singular values. */

        /*for (i = 0; i < n; i++) R.Vt.value[x]xv1[id++] = w1[i];*/
        /*id += nrow;*/
        nsig++;
      }
      id2++;
    }

    /* Rotate the singular vectors and values. */
    /* x is now the location of the highest singular value. */
    rotateArray(R.Vt.value[0], R.Vt.rows*R.Vt.cols,
      x*R.Vt.cols);
    R.d = Math.min(R.d, nsig);
    for (x = 0; x < R.d; x++) {
      /* multiply by matrix B first */
      svd_opb(A, R.Vt.value[x], xv2, OPBTemp);
      tmp0 = svd_ddot(n, R.Vt.value[x], 1, xv2, 1);
      svd_daxpy(n, -tmp0, R.Vt.value[x], 1, xv2, 1);
      tmp0 = sqrt(tmp0);
      xnorm = sqrt(svd_ddot(n, xv2, 1, xv2, 1));

      /* multiply by matrix A to get (scaled) left s-vector */
      svd_opa(A, R.Vt.value[x], R.Ut.value[x]);
      tmp1 = 1.0/tmp0;
      svd_dscal(A.rows, tmp1, R.Ut.value[x], 1);
      xnorm *= tmp1;
      bnd[i] = xnorm;
      R.S[x] = tmp0;
    }
  }
  SAFE_FREE(s);
  SAFE_FREE(xv2);
  return nsig;
}

/***********************************************************************
 *                                                                     *
 *                          lanso()                                    *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

 Description
 -----------

 Function determines when the restart of the Lanczos algorithm should
 occur and when it should terminate.

 Arguments
 ---------

 (input)
 n         dimension of the eigenproblem for matrix B
 iterations    upper limit of desired number of lanczos steps
 dimensions    upper limit of desired number of eigenpairs
 endl      left end of interval containing unwanted eigenvalues
 endr      right end of interval containing unwanted eigenvalues
 ritz      array to hold the ritz values
 bnd       array to hold the error bounds
 wptr      array of pointers that point to work space:
 wptr[0]-wptr[5]  six vectors of length n
 wptr[6] array to hold diagonal of the tridiagonal matrix T
 wptr[9] array to hold off-diagonal of T
 wptr[7] orthogonality estimate of Lanczos vectors at
 step j
 wptr[8] orthogonality estimate of Lanczos vectors at
 step j-1

 (output)
 j         number of Lanczos steps actually taken
 neig      number of ritz values stabilized
 ritz      array to hold the ritz values
 bnd       array to hold the error bounds
 ierr      (globally declared) error flag
 ierr = 8192 if stpone() fails to find a starting vector
 ierr = k if convergence did not occur for k-th eigenvalue
 in imtqlb()
 ierr = 0 otherwise


 Functions used
 --------------

 LAS    stpone, error_bound, lanczos_step
 MISC    svd_dsort2
 UTILITY  Math.min, svd_imax

 ***********************************************************************/

//int lanso(SMat A, long iterations, long dimensions, double endl,
//          double endr, double *ritz, double *bnd, double *wptr[],
//          long n) {
export function lanso(A, iterations, dimensions, endl, endr, ritz, bnd, wptr, n) {
  let alf, eta, oldeta, bet, wrk;
  let ll, first, last, ENOUGH, id2, id3, i, l, neig, j = 0, intro = 0;
  let neigp;

  alf = wptr[6];
  eta = wptr[7];
  oldeta = wptr[8];
  bet = wptr[9];
  wrk = wptr[5];

  /* take the first step */
  let {rnm, tol} = stpone(A, wptr, n);

  if (!rnm || ierr) return 0;
  eta[0] = eps1;
  oldeta[0] = eps1;
  ll = 0;
  first = 1;
  last = Math.min(dimensions + svd_imax(8, dimensions), iterations);
  ENOUGH = false;
  /*id1 = 0;*/
  while (/*id1 < dimensions && */!ENOUGH) {
    if (rnm <= tol) rnm = 0.0;

    /* the actual lanczos loop */
    let ret2 = lanczos_step(A, first, last, wptr, alf, eta, oldeta, bet, ll, ENOUGH, rnm, tol, n);

    j = ret2.j;
    rnm = ret2.rnm;
    tol = ret2.tol;
    ENOUGH = ret2.ENOUGH;
    ll = ret2.ll;

    if (ENOUGH) j = j - 1;
    else j = last - 1;
    first = j + 1;
    bet[j + 1] = rnm;

    /* analyze T */
    l = 0;
    for (id2 = 0; id2 < j; id2++) {
      if (l > j) break;
      for (i = l; i <= j; i++) {
        if (!bet[i + 1]) break;
      }
      if (i > j) i = j;

      /* now i is at the end of an unreduced submatrix */
      svd_dcopy(i - l + 1, alf, 1, ritz, -1, l, l);
      svd_dcopy(i - l, bet, 1, wrk, -1, l + 1, l + 1);

      imtqlb(i - l + 1, ritz, wrk, bnd, l);

      if (ierr) {
        printf("svdLAS2: imtqlb failed to converge (ierr = %ld)\n", ierr);
        printf("  l = %ld  i = %ld\n", l, i);

        for (id3 = l; id3 <= i; id3++) {
          printf("  %ld  %lg  %lg  %lg\n",
            id3, ritz[id3], wrk[id3], bnd[id3]);
        }
      }
      for (id3 = l; id3 <= i; id3++) {
        bnd[id3] = rnm*fabs(bnd[id3]);
      }
      l = i + 1;
    }

    /* sort eigenvalues into increasing order */
    svd_dsort2((j + 1)/2, j + 1, ritz, bnd);

    /*    for (i = 0; i < iterations; i++)
      printf("%f ", ritz[i]);
      printf("\n"); */

    /* massage error bounds for very close ritz values */
    let ret = error_bound(ENOUGH, endl, endr, ritz, bnd, j, tol);

    neig = ret.neig;
    ENOUGH = ret.enough;

    /* should we stop? */
    if (neig < dimensions) {
      if (!neig) {
        last = first + 9;
        intro = first;
      } else last = first + svd_imax(3, 1 + ((j - intro)*(dimensions - neig))/
        neig);
      last = Math.min(last, iterations);
    } else ENOUGH = true;
    ENOUGH = ENOUGH || first >= iterations;
    /* id1++; */
    /* printf("id1=%d dimen=%d first=%d\n", id1, dimensions, first); */
  }
  store(n, STORQ, j, wptr[1]);
  return {steps: j, neig};
}


/***********************************************************************
 *                                                                     *
 *      lanczos_step()                                 *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

 Description
 -----------

 Function embodies a single Lanczos step

 Arguments
 ---------

 (input)
 n        dimension of the eigenproblem for matrix B
 first    start of index through loop
 last     end of index through loop
 wptr      array of pointers pointing to work space
 alf      array to hold diagonal of the tridiagonal matrix T
 eta      orthogonality estimate of Lanczos vectors at step j
 oldeta   orthogonality estimate of Lanczos vectors at step j-1
 bet      array to hold off-diagonal of T
 ll       number of intitial Lanczos vectors in local orthog.
 (has value of 0, 1 or 2)
 enough   stop flag

 Functions used
 --------------

 BLAS    svd_ddot, svd_dscal, svd_daxpy, svd_datx, svd_dcopy
 USER    store
 LAS    purge, ortbnd, startv
 UTILITY  Math.min, svd_imax

 ***********************************************************************/

/*long lanczos_step(SMat A, long first, long last, double *wptr[],
		  double *alf, double *eta, double *oldeta,
		  double *bet, long *ll, long *enough, double *rnmp,
                  double *tolp, long n) {*/
export function lanczos_step(A, first, last, wptr,
                             alf, eta, oldeta, bet, ll, enough, rnmp, tolp, n) {

  let t, mid, rnm = rnmp, tol = tolp, anorm;
  let i, j;

  for (j = first; j < last; j++) {
    mid = wptr[2];
    wptr[2] = wptr[1];
    wptr[1] = mid;
    mid = wptr[3];
    wptr[3] = wptr[4];
    wptr[4] = mid;

    store(n, STORQ, j - 1, wptr[2]);
    if (j - 1 < MAXLL) store(n, STORP, j - 1, wptr[4]);
    bet[j] = rnm;

    /* restart if invariant subspace is found */
    if (!bet[j]) {
      rnm = startv(A, wptr, j, n);
      if (ierr) return j;
      if (!rnm)
        enough = true;
    }
    if (enough) {
      /* added by Doug... */
      /* These lines fix a bug that occurs with low-rank matrices */
      mid = wptr[2];
      wptr[2] = wptr[1];
      wptr[1] = mid;
      /* ...added by Doug */
      break;
    }

    /* take a lanczos step */
    t = 1.0/rnm;
    svd_datx(n, t, wptr[0], 1, wptr[1], 1);
    svd_dscal(n, t, wptr[3], 1);
    svd_opb(A, wptr[3], wptr[0], OPBTemp);
    svd_daxpy(n, -rnm, wptr[2], 1, wptr[0], 1);
    alf[j] = svd_ddot(n, wptr[0], 1, wptr[3], 1);
    svd_daxpy(n, -alf[j], wptr[1], 1, wptr[0], 1);

    /* orthogonalize against initial lanczos vectors */
    if (j <= MAXLL && (fabs(alf[j - 1]) > 4.0*fabs(alf[j]))) {
      ll = j;
    }

    for (i = 0; i < Math.min(ll, j - 1); i++) {
      store(n, RETRP, i, wptr[5]);
      t = svd_ddot(n, wptr[5], 1, wptr[0], 1);
      store(n, RETRQ, i, wptr[5]);
      svd_daxpy(n, -t, wptr[5], 1, wptr[0], 1);
      eta[i] = eps1;
      oldeta[i] = eps1;
    }

    /* extended local reorthogonalization */
    t = svd_ddot(n, wptr[0], 1, wptr[4], 1);
    svd_daxpy(n, -t, wptr[2], 1, wptr[0], 1);
    if (bet[j] > 0.0)
      bet[j] = bet[j] + t;
    t = svd_ddot(n, wptr[0], 1, wptr[3], 1);
    svd_daxpy(n, -t, wptr[1], 1, wptr[0], 1);
    alf[j] = alf[j] + t;
    svd_dcopy(n, wptr[0], 1, wptr[4], 1);
    rnm = sqrt(svd_ddot(n, wptr[0], 1, wptr[4], 1));
    anorm = bet[j] + fabs(alf[j]) + rnm;
    tol = reps*anorm;

    /* update the orthogonality bounds */
    ortbnd(alf, eta, oldeta, bet, j, rnm);

    /* restore the orthogonality state when needed */
    rnm = purge(n, ll, wptr[0], wptr[1], wptr[4], wptr[3], wptr[5], eta, oldeta,
      j, rnm, tol);
    if (rnm <= tol) rnm = 0.0;
  }

  return {j, rnm, tol, ENOUGH: enough, ll};
}

/***********************************************************************
 *                                                                     *
 *                          ortbnd()                                   *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

 Description
 -----------

 Funtion updates the eta recurrence

 Arguments
 ---------

 (input)
 alf      array to hold diagonal of the tridiagonal matrix T
 eta      orthogonality estimate of Lanczos vectors at step j
 oldeta   orthogonality estimate of Lanczos vectors at step j-1
 bet      array to hold off-diagonal of T
 n        dimension of the eigenproblem for matrix B
 j        dimension of T
 rnm      norm of the next residual vector
 eps1      roundoff estimate for dot product of two unit vectors

 (output)
 eta      orthogonality estimate of Lanczos vectors at step j+1
 oldeta   orthogonality estimate of Lanczos vectors at step j


 Functions used
 --------------

 BLAS    svd_dswap

 ***********************************************************************/

//void ortbnd(double *alf, double *eta, double *oldeta, double *bet, long step,
//            double rnm) {
export function ortbnd(alf, eta, oldeta, bet, step, rnm) {
  let i;
  if (step < 1) return;
  if (rnm) {
    if (step > 1) {
      oldeta[0] = (bet[1]*eta[1] + (alf[0] - alf[step])*eta[0] -
        bet[step]*oldeta[0])/rnm + eps1;
    }
    for (i = 1; i <= step - 2; i++) {
      oldeta[i] = (bet[i + 1]*eta[i + 1] + (alf[i] - alf[step])*eta[i] +
        bet[i]*eta[i - 1] - bet[step]*oldeta[i])/rnm + eps1;
    }
  }

  oldeta[step - 1] = eps1;
  svd_dswap(step, oldeta, 1, eta, 1);
  eta[step] = eps1;
  return;
}

/***********************************************************************
 *                                                                     *
 *        purge()                                *
 *                                                                     *
 ***********************************************************************/

/***********************************************************************

 Description
 -----------

 Function examines the state of orthogonality between the new Lanczos
 vector and the previous ones to decide whether re-orthogonalization
 should be performed


 Arguments
 ---------

 (input)
 n        dimension of the eigenproblem for matrix B
 ll       number of intitial Lanczos vectors in local orthog.
 r        residual vector to become next Lanczos vector
 q        current Lanczos vector
 ra       previous Lanczos vector
 qa       previous Lanczos vector
 wrk      temporary vector to hold the previous Lanczos vector
 eta      state of orthogonality between r and prev. Lanczos vectors
 oldeta   state of orthogonality between q and prev. Lanczos vectors
 j        current Lanczos step

 (output)
 r      residual vector orthogonalized against previous Lanczos
 vectors
 q        current Lanczos vector orthogonalized against previous ones


 Functions used
 --------------

 BLAS    svd_daxpy,  svd_dcopy,  svd_idamax,  svd_ddot
 USER    store

 ***********************************************************************/

/*****************************************************************
 * Function finds the index of element having max. absolute value*
 * based on FORTRAN 77 routine from Linpack by J. Dongarra       *
 *****************************************************************/
//long svd_idamax(long n, double *dx, long incx) {
export function svd_idamax(n, dx, incx, start = 0) {
  let ix, i, imax, dtemp, dmax;

  if (n < 1) return (-1);
  if (n === 1) return (0);
  if (incx === 0) return (-1);

  if (incx < 0) ix = (-n + 1)*incx;
  else ix = 0;

  let di = start;

  imax = ix;
  di += ix;
  dmax = fabs(dx[di]);

  for (i = 1; i < n; i++) {
    ix += incx;
    di += incx;
    dtemp = fabs(dx[di]);
    if (dtemp > dmax) {
      dmax = dtemp;
      imax = ix;
    }
  }
  return (imax);
}

/*void purge(long n, long ll, double *r, double *q, double *ra,
     double *qa, double *wrk, double *eta, double *oldeta, long step,
           double *rnmp, double tol) {
*/
export function purge(n, ll, r, q, ra, qa, wrk, eta, oldeta, step, rnmp, tol) {
  let t, tq, tr, reps1, rnm = rnmp;
  let k, iteration, flag, i;

  if (step < ll + 2) {
    return;
  }

  k = svd_idamax(step - (ll + 1), eta, 1, ll) + ll;

  if (fabs(eta[k]) > reps) {
    reps1 = eps1/reps;
    iteration = 0;
    flag = true;
    while (iteration < 2 && flag) {
      if (rnm > tol) {

        /* bring in a lanczos vector t and orthogonalize both
         * r and q against it */
        tq = 0.0;
        tr = 0.0;
        for (i = ll; i < step; i++) {
          store(n, RETRQ, i, wrk);
          t = -svd_ddot(n, qa, 1, wrk, 1);
          tq += fabs(t);
          svd_daxpy(n, t, wrk, 1, q, 1);
          t = -svd_ddot(n, ra, 1, wrk, 1);
          tr += fabs(t);
          svd_daxpy(n, t, wrk, 1, r, 1);
        }
        svd_dcopy(n, q, 1, qa, 1);
        t = -svd_ddot(n, r, 1, qa, 1);
        tr += fabs(t);
        svd_daxpy(n, t, q, 1, r, 1);
        svd_dcopy(n, r, 1, ra, 1);
        rnm = sqrt(svd_ddot(n, ra, 1, r, 1));
        if (tq <= reps1 && tr <= reps1*rnm) flag = false;
      }
      iteration++;
    }
    for (i = ll; i <= step; i++) {
      eta[i] = eps1;
      oldeta[i] = eps1;
    }
  }

  return rnm;
}


/***********************************************************************
 *                                                                     *
 *                         stpone()                                    *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

 Description
 -----------

 Function performs the first step of the Lanczos algorithm.  It also
 does a step of extended local re-orthogonalization.

 Arguments
 ---------

 (input)
 n      dimension of the eigenproblem for matrix B

 (output)
 ierr   error flag
 wptr   array of pointers that point to work space that contains
 wptr[0]             r[j]
 wptr[1]             q[j]
 wptr[2]             q[j-1]
 wptr[3]             p
 wptr[4]             p[j-1]
 wptr[6]             diagonal elements of matrix T


 Functions used
 --------------

 BLAS    svd_daxpy, svd_datx, svd_dcopy, svd_ddot, svd_dscal
 USER    store, opb
 LAS    startv

 ***********************************************************************/

//void stpone(SMat A, double *wrkptr[], double *rnmp, double *tolp, long n) {
export function stpone(A, wrkptr, n) {
  let t, alf, rnm, anorm;
  alf = wrkptr[6];

  /* get initial vector; default is random */
  rnm = startv(A, wrkptr, 0, n);
  if (rnm === 0.0 || ierr !== 0) return;

  /* normalize starting vector */
  t = 1.0/rnm;
  svd_datx(n, t, wrkptr[0], 1, wrkptr[1], 1);
  svd_dscal(n, t, wrkptr[3], 1);

  /* take the first step */
  svd_opb(A, wrkptr[3], wrkptr[0], OPBTemp);
  alf[0] = svd_ddot(n, wrkptr[0], 1, wrkptr[3], 1);
  svd_daxpy(n, -alf[0], wrkptr[1], 1, wrkptr[0], 1);
  t = svd_ddot(n, wrkptr[0], 1, wrkptr[3], 1);
  svd_daxpy(n, -t, wrkptr[1], 1, wrkptr[0], 1);
  alf[0] += t;
  svd_dcopy(n, wrkptr[0], 1, wrkptr[4], 1);
  rnm = sqrt(svd_ddot(n, wrkptr[0], 1, wrkptr[4], 1));
  anorm = rnm + fabs(alf[0]);

  let tol = reps*anorm;

  return {rnm, tol};
}

/***********************************************************************
 *                                                                     *
 *                         startv()                                    *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

 Description
 -----------

 Function delivers a starting vector in r and returns |r|; it returns
 zero if the range is spanned, and ierr is non-zero if no starting
 vector within range of operator can be found.

 Parameters
 ---------

 (input)
 n      dimension of the eigenproblem matrix B
 wptr   array of pointers that point to work space
 j      starting index for a Lanczos run
 eps    machine epsilon (relative precision)

 (output)
 wptr   array of pointers that point to work space that contains
 r[j], q[j], q[j-1], p[j], p[j-1]
 ierr   error flag (nonzero if no starting vector can be found)

 Functions used
 --------------

 BLAS    svd_ddot, svd_dcopy, svd_daxpy
 USER    svd_opb, store
 MISC    random

 ***********************************************************************/

export class RandState {
  constructor() {
    this.i = 0;
  }
}

let randstates = new cachering(() => new RandState(), 4196);

//double startv(SMat A, double *wptr[], long step, long n) {
function startv(A, wptr, step, n) {
  let rnm2, r, t;
  let irand = randstates.next();
  let id, i;

  /* get initial vector; default is random */
  rnm2 = svd_ddot(n, wptr[0], 1, wptr[0], 1);
  irand.i = 918273 + step;
  r = wptr[0];
  for (id = 0; id < 3; id++) {
    if (id > 0 || step > 0 || rnm2 === 0)
      for (i = 0; i < n; i++) {
        r[i] = svd_random2(irand);
      }
    svd_dcopy(n, wptr[0], 1, wptr[3], 1);

    /* apply operator to put r in range (essential if m singular) */
    svd_opb(A, wptr[3], wptr[0], OPBTemp);
    svd_dcopy(n, wptr[0], 1, wptr[3], 1);
    rnm2 = svd_ddot(n, wptr[0], 1, wptr[3], 1);
    if (rnm2 > 0.0) break;
  }

  /* fatal error */
  if (rnm2 <= 0.0) {
    ierr = 8192;
    return (-1);
  }
  if (step > 0) {
    for (i = 0; i < step; i++) {
      store(n, RETRQ, i, wptr[5]);
      t = -svd_ddot(n, wptr[3], 1, wptr[5], 1);
      svd_daxpy(n, t, wptr[5], 1, wptr[0], 1);
    }

    /* make sure q[step] is orthogonal to q[step-1] */
    t = svd_ddot(n, wptr[4], 1, wptr[0], 1);
    svd_daxpy(n, -t, wptr[2], 1, wptr[0], 1);
    svd_dcopy(n, wptr[0], 1, wptr[3], 1);
    t = svd_ddot(n, wptr[3], 1, wptr[0], 1);
    if (t <= eps*rnm2) t = 0.0;
    rnm2 = t;
  }
  return (sqrt(rnm2));
}

/***********************************************************************
 *                                                                     *
 *      error_bound()                                  *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

 Description
 -----------

 Function massages error bounds for very close ritz values by placing
 a gap between them.  The error bounds are then refined to reflect
 this.


 Arguments
 ---------

 (input)
 endl     left end of interval containing unwanted eigenvalues
 endr     right end of interval containing unwanted eigenvalues
 ritz     array to store the ritz values
 bnd      array to store the error bounds
 enough   stop flag


 Functions used
 --------------

 BLAS    svd_idamax
 UTILITY  svd_dmin

 ***********************************************************************/

//long error_bound(long *enough, double endl, double endr,
//                 double *ritz, double *bnd, long step, double tol) {
export function error_bound(enoughIn, endl, endr, ritz, bnd, step, tol) {
  let mid, i, neig, enough = enoughIn;
  let gapl, gap;

  /* massage error bounds for very close ritz values */
  mid = svd_idamax(step + 1, bnd, 1);

  for (i = ((step + 1) + (step - 1))/2; i >= mid + 1; i -= 1) {
    if (fabs(ritz[i - 1] - ritz[i]) < eps34*fabs(ritz[i]))
      if (bnd[i] > tol && bnd[i - 1] > tol) {
        bnd[i - 1] = sqrt(bnd[i]*bnd[i] + bnd[i - 1]*bnd[i - 1]);
        bnd[i] = 0.0;
      }
  }


  for (i = ((step + 1) - (step - 1))/2; i <= mid - 1; i += 1) {
    if (fabs(ritz[i + 1] - ritz[i]) < eps34*fabs(ritz[i]))
      if (bnd[i] > tol && bnd[i + 1] > tol) {
        bnd[i + 1] = sqrt(bnd[i]*bnd[i] + bnd[i + 1]*bnd[i + 1]);
        bnd[i] = 0.0;
      }
  }

  /* refine the error bounds */
  neig = 0;
  gapl = ritz[step] - ritz[0];
  for (i = 0; i <= step; i++) {
    gap = gapl;
    if (i < step) gapl = ritz[i + 1] - ritz[i];
    gap = svd_dmin(gap, gapl);
    if (gap > bnd[i]) bnd[i] = bnd[i]*(bnd[i]/gap);
    if (bnd[i] <= 16.0*eps*fabs(ritz[i])) {
      neig++;
      if (!enough) enough = endl < ritz[i] && ritz[i] < endr;
    }
  }
  return {neig, enough};
}

/***********************************************************************
 *                                                                     *
 *        imtqlb()             *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

 Description
 -----------

 imtqlb() is a translation of a Fortran version of the Algol
 procedure IMTQL1, Num. Math. 12, 377-383(1968) by Martin and
 Wilkinson, as modified in Num. Math. 15, 450(1970) by Dubrulle.
 Handbook for Auto. Comp., vol.II-Linear Algebra, 241-248(1971).
 See also B. T. Smith et al, Eispack Guide, Lecture Notes in
 Computer Science, Springer-Verlag, (1976).

 The function finds the eigenvalues of a symmetric tridiagonal
 matrix by the implicit QL method.


 Arguments
 ---------

 (input)
 n      order of the symmetric tridiagonal matrix
 d      contains the diagonal elements of the input matrix
 e      contains the subdiagonal elements of the input matrix in its
 last n-1 positions.  e[0] is arbitrary

 (output)
 d      contains the eigenvalues in ascending order.  if an error
 exit is made, the eigenvalues are correct and ordered for
 indices 0,1,...ierr, but may not be the smallest eigenvalues.
 e      has been destroyed.
 ierr   set to zero for normal return, j if the j-th eigenvalue has
 not been determined after 30 iterations.

 Functions used
 --------------

 UTILITY  svd_fsign
 MISC    svd_pythag

 ***********************************************************************/

//void imtqlb(long n, double d[], double e[], double bnd[])
export function imtqlb(n, d, e, bnd, start = 0) {
  let last, l, m, i, iteration;

  /* various flags */
  let exchange, convergence, underflow;

  let b, test, g, r, s, c, p, f;

  if (n === 1) return;
  ierr = 0;
  bnd[start] = 1.0;
  last = n - 1;
  for (i = 1; i < n; i++) {
    bnd[i + start] = 0.0;
    e[i - 1 + start] = e[i + start];
  }
  e[last + start] = 0.0;
  for (l = 0; l < n; l++) {
    iteration = 0;
    while (iteration <= 30) {
      for (m = l; m < n; m++) {
        convergence = false;
        if (m === last) break;
        else {
          test = fabs(d[m + start]) + fabs(d[m + 1 + start]);
          if (test + fabs(e[m + start]) === test) convergence = true;
        }
        if (convergence) break;
      }
      p = d[l + start];
      f = bnd[l + start];
      if (m !== l) {
        if (iteration === 30) {
          ierr = l;
          return;
        }
        iteration += 1;
        /*........ form shift ........*/
        g = (d[l + 1 + start] - p)/(2.0*e[l + start]);
        r = svd_pythag(g, 1.0);
        g = d[m + start] - p + e[l + start]/(g + svd_fsign(r, g));
        s = 1.0;
        c = 1.0;
        p = 0.0;
        underflow = false;
        i = m - 1;
        while (underflow === false && i >= l) {
          f = s*e[i + start];
          b = c*e[i + start];
          r = svd_pythag(f, g);
          e[i + 1 + start] = r;
          if (r === 0.0) underflow = true;
          else {
            s = f/r;
            c = g/r;
            g = d[i + 1 + start] - p;
            r = (d[i + start] - g)*s + 2.0*c*b;
            p = s*r;
            d[i + 1 + start] = g + p;
            g = c*r - b;
            f = bnd[i + 1 + start];
            bnd[i + 1 + start] = s*bnd[i + start] + c*f;
            bnd[i + start] = c*bnd[i + start] - s*f;
            i--;
          }
        }       /* end while (underflow !== false && i >= l) */
        /*........ recover from underflow .........*/
        if (underflow) {
          d[i + 1 + start] -= p;
          e[m + start] = 0.0;
        } else {
          d[l + start] -= p;
          e[l + start] = g;
          e[m + start] = 0.0;
        }
      } 		       		   /* end if (m !== l) */
      else {

        /* order the eigenvalues */
        exchange = true;
        if (l !== 0) {
          i = l;
          while (i >= 1 && exchange === true) {
            if (p < d[i - 1 + start]) {
              d[i + start] = d[i - 1 + start];
              bnd[i + start] = bnd[i - 1 + start];
              i--;
            } else exchange = false;
          }
        }
        if (exchange) i = 0;
        d[i + start] = p;
        bnd[i + start] = f;
        iteration = 31;
      }
    }			       /* end while (iteration <= 30) */
  }				   /* end for (l=0; l<n; l++) */
  return;
}						  /* end main */

/***********************************************************************
 *                                                                     *
 *        imtql2()             *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

 Description
 -----------

 imtql2() is a translation of a Fortran version of the Algol
 procedure IMTQL2, Num. Math. 12, 377-383(1968) by Martin and
 Wilkinson, as modified in Num. Math. 15, 450(1970) by Dubrulle.
 Handbook for Auto. Comp., vol.II-Linear Algebra, 241-248(1971).
 See also B. T. Smith et al, Eispack Guide, Lecture Notes in
 Computer Science, Springer-Verlag, (1976).

 This function finds the eigenvalues and eigenvectors of a symmetric
 tridiagonal matrix by the implicit QL method.


 Arguments
 ---------

 (input)
 nm     row dimension of the symmetric tridiagonal matrix
 n      order of the matrix
 d      contains the diagonal elements of the input matrix
 e      contains the subdiagonal elements of the input matrix in its
 last n-1 positions.  e[0] is arbitrary
 z      contains the identity matrix

 (output)
 d      contains the eigenvalues in ascending order.  if an error
 exit is made, the eigenvalues are correct but unordered for
 for indices 0,1,...,ierr.
 e      has been destroyed.
 z      contains orthonormal eigenvectors of the symmetric
 tridiagonal (or full) matrix.  if an error exit is made,
 z contains the eigenvectors associated with the stored
 eigenvalues.
 ierr   set to zero for normal return, j if the j-th eigenvalue has
 not been determined after 30 iterations.


 Functions used
 --------------
 UTILITY  svd_fsign
 MISC    svd_pythag

 ***********************************************************************/

//void imtql2(long nm, long n, double d[], double e[], double z[])
export function imtql2(nm, n, d, e, z) {
  let index, nnm, j, last, l, m, i, k, iteration, convergence, underflow;
  let b, test, g, r, s, c, p, f;

  if (n === 1) return;
  ierr = 0;
  last = n - 1;
  for (i = 1; i < n; i++) {
    e[i - 1] = e[i];
  }
  e[last] = 0.0;
  nnm = n*nm;
  for (l = 0; l < n; l++) {
    iteration = 0;

    /* look for small sub-diagonal element */
    while (iteration <= 30) {
      for (m = l; m < n; m++) {
        convergence = false;
        if (m === last) break;
        else {
          test = fabs(d[m]) + fabs(d[m + 1]);
          if (test + fabs(e[m]) === test) convergence = true;
        }
        if (convergence) break;
      }
      if (m !== l) {

        /* set error -- no convergence to an eigenvalue after
         * 30 iterations. */
        if (iteration === 30) {
          ierr = l;
          return;
        }
        p = d[l];
        iteration += 1;

        /* form shift */
        g = (d[l + 1] - p)/(2.0*e[l]);
        r = svd_pythag(g, 1.0);
        g = d[m] - p + e[l]/(g + svd_fsign(r, g));
        s = 1.0;
        c = 1.0;
        p = 0.0;
        underflow = false;
        i = m - 1;
        while (underflow === false && i >= l) {
          f = s*e[i];
          b = c*e[i];
          r = svd_pythag(f, g);
          e[i + 1] = r;
          if (r === 0.0) underflow = true;
          else {
            s = f/r;
            c = g/r;
            g = d[i + 1] - p;
            r = (d[i] - g)*s + 2.0*c*b;
            p = s*r;
            d[i + 1] = g + p;
            g = c*r - b;

            /* form vector */
            for (k = 0; k < nnm; k += n) {
              index = k + i;
              f = z[index + 1];
              z[index + 1] = s*z[index] + c*f;
              z[index] = c*z[index] - s*f;
            }
            i--;
          }
        }   /* end while (underflow !== false && i >= l) */
        /*........ recover from underflow .........*/
        if (underflow) {
          d[i + 1] -= p;
          e[m] = 0.0;
        } else {
          d[l] -= p;
          e[l] = g;
          e[m] = 0.0;
        }
      } else break;
    }		/*...... end while (iteration <= 30) .........*/
  }		/*...... end for (l=0; l<n; l++) .............*/

  /* order the eigenvalues */
  for (l = 1; l < n; l++) {
    i = l - 1;
    k = i;
    p = d[i];
    for (j = l; j < n; j++) {
      if (d[j] < p) {
        k = j;
        p = d[j];
      }
    }
    /* ...and corresponding eigenvectors */
    if (k !== i) {
      d[k] = d[i];
      d[i] = p;
      for (j = 0; j < nnm; j += n) {
        p = z[j + i];
        z[j + i] = z[j + k];
        z[j + k] = p;
      }
    }
  }
  return;
}		/*...... end main ............................*/

/***********************************************************************
 *                                                                     *
 *        machar()             *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

 Description
 -----------

 This function is a partial translation of a Fortran-77 subroutine
 written by W. J. Cody of Argonne National Laboratory.
 It dynamically determines the listed machine parameters of the
 floating-point arithmetic.  According to the documentation of
 the Fortran code, "the determination of the first three uses an
 extension of an algorithm due to M. Malcolm, ACM 15 (1972),
 pp. 949-951, incorporating some, but not all, of the improvements
 suggested by M. Gentleman and S. Marovich, CACM 17 (1974),
 pp. 276-277."  The complete Fortran version of this translation is
 documented in W. J. Cody, "Machar: a Subroutine to Dynamically
 Determine Determine Machine Parameters," TOMS 14, December, 1988.


 Parameters reported
 -------------------

 ibeta     the radix for the floating-point representation
 it        the number of base ibeta digits in the floating-point
 significand
 irnd      0 if floating-point addition chops
 1 if floating-point addition rounds, but not in the
 ieee style
 2 if floating-point addition rounds in the ieee style
 3 if floating-point addition chops, and there is
 partial underflow
 4 if floating-point addition rounds, but not in the
 ieee style, and there is partial underflow
 5 if floating-point addition rounds in the ieee style,
 and there is partial underflow
 machep    the largest negative integer such that
 1.0+float(ibeta)**machep .ne. 1.0, except that
 machep is bounded below by  -(it+3)
 negeps    the largest negative integer such that
 1.0-float(ibeta)**negeps .ne. 1.0, except that
 negeps is bounded below by  -(it+3)

 ***********************************************************************/

let machar_rets = new cachering(() => {
  return {
    ibeta: 0, it: 0, irnd: 0, machep: 0, negep: 0, eps0: 0
  }
}, 64);

export function machar() {
//void machar(long *ibeta, long *it, long *irnd, long *machep, long *negep) {
  let ret = machar_rets.next();

  let beta, betain, betah, a, b, ZERO, ONE, TWO, temp, tempa, temp1;
  let i, itemp;

  ONE = 1.0;
  TWO = ONE + ONE;
  ZERO = ONE - ONE;

  a = ONE;
  temp1 = ONE;
  while (temp1 - ONE === ZERO) {
    a = a + a;
    temp = a + ONE;
    temp1 = temp - a;
    b += a; /* to prevent icc compiler error */
  }
  b = ONE;
  itemp = 0;
  while (itemp === 0) {
    b = b + b;
    temp = a + b;
    itemp = (long)(temp - a);
  }

  ret.ibeta = itemp;
  beta = ret.ibeta;

  ret.it = 0;
  b = ONE;
  temp1 = ONE;
  while (temp1 - ONE === ZERO) {
    ret.it = ret.it + 1;
    b = b*beta;
    temp = b + ONE;
    temp1 = temp - b;
  }
  ret.irnd = 0;
  betah = beta/TWO;
  temp = a + betah;
  if (temp - a !== ZERO) ret.irnd = 1;
  tempa = a + beta;
  temp = tempa + betah;
  if ((ret.irnd === 0) && (temp - tempa !== ZERO)) ret.irnd = 2;

  ret.negep = ret.it + 3;
  betain = ONE/beta;
  a = ONE;
  for (i = 0; i < ret.negep; i++) {
    a = a*betain;
  }
  b = a;
  temp = ONE - a;
  while (temp - ONE === ZERO) {
    a = a*beta;
    ret.negep = ret.negep - 1;
    temp = ONE - a;
  }
  ret.negep = -(ret.negep);

  ret.machep = -(ret.it) - 3;
  a = b;
  temp = ONE + a;
  while (temp - ONE === ZERO) {
    a = a*beta;
    ret.machep = ret.machep + 1;
    temp = ONE + a;
  }

  ret.eps0 = a;

  return ret;
}

/***********************************************************************
 *                                                                     *
 *                     store()                                         *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

 Description
 -----------

 store() is a user-supplied function which, based on the input
 operation flag, stores to or retrieves from memory a vector.


 Arguments
 ---------

 (input)
 n       length of vector to be stored or retrieved
 isw     operation flag:
 isw = 1 request to store j-th Lanczos vector q(j)
 isw = 2 request to retrieve j-th Lanczos vector q(j)
 isw = 3 request to store q(j) for j = 0 or 1
 isw = 4 request to retrieve q(j) for j = 0 or 1
 s     contains the vector to be stored for a "store" request

 (output)
 s     contains the vector retrieved for a "retrieve" request

 Functions used
 --------------

 BLAS    svd_dcopy

 ***********************************************************************/

export function store(n, isw, j, s) {
  /* printf("called store %ld %ld\n", isw, j); */
  switch (isw) {
    case STORQ:
      if (!LanStore[j + MAXLL]) {
        if (!(LanStore[j + MAXLL] = svd_doubleArray(n, false, "LanStore[j]")))
          throw new SvdError("svdLAS2: failed to allocate LanStore[%d]", j + MAXLL);
      }
      svd_dcopy(n, s, 1, LanStore[j + MAXLL], 1);
      break;
    case RETRQ:
      if (!LanStore[j + MAXLL])
        throw new SvdError("svdLAS2: store (RETRQ) called on index %d (not allocated)",
          j + MAXLL);
      svd_dcopy(n, LanStore[j + MAXLL], 1, s, 1);
      break;
    case STORP:
      if (j >= MAXLL) {
        throw new SvdError("svdLAS2: store (STORP) called with j >= MAXLL");
        break;
      }
      if (!LanStore[j]) {
        if (!(LanStore[j] = svd_doubleArray(n, false, "LanStore[j]")))
          throw new SvdError("svdLAS2: failed to allocate LanStore[%d]", j);
      }
      svd_dcopy(n, s, 1, LanStore[j], 1);
      break;
    case RETRP:
      if (j >= MAXLL) {
        throw new SvdError("svdLAS2: store (RETRP) called with j >= MAXLL");
        break;
      }
      if (!LanStore[j])
        throw new SvdError("svdLAS2: store (RETRP) called on index %d (not allocated)",
          j);
      svd_dcopy(n, LanStore[j], 1, s, 1);
      break;
  }
}


/* Row major order.  Rows are vectors that are consecutive in memory.  Matrix
   is initialized to empty. */
export function svdNewDMat(rows, cols) {
  let i;
  let D = new DMat();

  D.rows = rows;
  D.cols = cols;

  D.value = tempArray(rows, true);
  D.value[0] = dtempArray(rows*cols);

  for (i = 1; i < rows; i++) {
    D.value[i] = D.value[i - 1] + cols;
  }

  return D;
}

export function svdFreeDMat(D) {
}


export function svdNewSMat(rows, cols, vals) {
  let S = new SMat();

  S.rows = rows;
  S.cols = cols;
  S.vals = vals;
  S.pointr = svd_longArray(cols + 1, true, "svdNewSMat: pointr");
  S.rowind = svd_longArray(vals, false, "svdNewSMat: rowind");
  S.value = svd_doubleArray(vals, false, "svdNewSMat: value");
  return S;
}

export function svdFreeSMat(S) {
  if (!S) {
    return;
  }

  SAFE_FREE(S.pointr);
  SAFE_FREE(S.rowind);
  SAFE_FREE(S.value);

  //free(S);
}


/* Creates an empty SVD record */
export function svdNewSVDRec() {
  return new SVDRec();
}

/* Frees an svd rec and all its contents. */
export function svdFreeSVDRec(R) {
}


/**************************** Conversion *************************************/

/* Converts a sparse matrix to a dense one (without affecting the former) */
export function svdConvertStoD(S) {
  let i, c;
  let D = svdNewDMat(S.rows, S.cols);

  for (i = 0, c = 0; i < S.vals; i++) {
    while (S.pointr[c + 1] <= i) {
      c++;
    }
    D.value[S.rowind[i]][c] = S.value[i];
  }

  return D;
}

/* Converts a dense matrix to a sparse one (without affecting the dense one) */
export function svdConvertDtoS(D) {
  let S;

  let i, j, n;
  for (i = 0, n = 0; i < D.rows; i++) {
    for (j = 0; j < D.cols; j++) {
      if (D.value[i][j] != 0) n++;
    }
  }

  S = svdNewSMat(D.rows, D.cols, n);

  for (j = 0, n = 0; j < D.cols; j++) {
    S.pointr[j] = n;
    for (i = 0; i < D.rows; i++) {
      if (D.value[i][j] != 0) {
        S.rowind[n] = i;
        S.value[n] = D.value[i][j];
        n++;
      }
    }
  }
  S.pointr[S.cols] = S.vals;
  return S;
}

/* Transposes a dense matrix. */
export function svdTransposeD(D) {
  let r, c;
  let N = svdNewDMat(D.cols, D.rows);

  for (r = 0; r < D.rows; r++) {
    for (c = 0; c < D.cols; c++) {
      N.value[c][r] = D.value[r][c];
    }
  }
  return N;
}

/* Efficiently transposes a sparse matrix. */
export function svdTransposeS(S) {
  let r, c, i, j;
  let N = svdNewSMat(S.cols, S.rows, S.vals);

  /* Count number nz in each row. */
  for (i = 0; i < S.vals; i++) {
    N.pointr[S.rowind[i]]++;
  }
  /* Fill each cell with the starting point of the previous row. */
  N.pointr[S.rows] = S.vals - N.pointr[S.rows - 1];
  for (r = S.rows - 1; r > 0; r--) {
    N.pointr[r] = N.pointr[r + 1] - N.pointr[r - 1];
  }
  N.pointr[0] = 0;
  /* Assign the new columns and values. */
  for (c = 0, i = 0; c < S.cols; c++) {
    for (; i < S.pointr[c + 1]; i++) {
      r = S.rowind[i];
      j = N.pointr[r + 1]++;
      N.rowind[j] = c;
      N.value[j] = S.value[i];
    }
  }
  return N;
}

