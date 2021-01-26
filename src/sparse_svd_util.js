let abs   = Math.abs,
    fabs  = Math.abs,
    sqrt  = Math.sqrt,
    floor = Math.floor,
    ceil  = Math.ceil,
    atan  = Math.atan,
    pow   = Math.pow,
    log   = Math.log;

export class SvdError extends Error {

}

/* Harwell-Boeing sparse matrix. */
export class SMat {
  constructor(pointr, rowind, value, rows, cols, vals) {
    this.rows = rows;
    this.cols = cols;
    this.vals = vals;
    this.pointr = pointr; /* For each col (plus 1), index of first non-zero entry. */
    this.rowind = rowind; /* For each nz entry, the row index. */
    this.value = value;  /* For each nz entry, the value. */
  }
}

export class SVDRec {
  constructor(d, Ut, S, Vt) {
    this.d = d;      /* Dimensionality (rank) */
    this.Ut = Ut;    /* Transpose of left singular vectors. (d by m)
                      The vectors are the rows of Ut. */
    this.S = S;      /* Array of singular values. (length d) */
    this.Vt = Vt;    /* Transpose of right singular vectors. (d by n)
                      The vectors are the rows of Vt. */
  }
};

/* Row-major dense matrix.  Rows are consecutive vectors. */
export class DMat {
  constructor(rows, cols, value) {
    this.rows = rows;
    this.cols = cols;
    this.value = vaule; /* Accessed by [row][col]. Free value[0] and value to free.*/
  }
}

export class cachering extends Array {
  constructor(func, num) {
    super();
    this.cur = 0;

    for (let i = 0; i < num; i++) {
      this.push(func());
    }
  }

  next() {
    let cur = this[this.cur];
    this.cur = (this.cur + 1)%this.length;

    return cur;
  }
}


export class ArrayPool {
  constructor(cls, clearval = 0.0) {
    this.pools = new Map();
    this.cls = cls;
    this.clearval = clearval;
  }

  get(n, clear) {
    if (!this.pools.has(n)) {
      this.pools.set(n, new cachering(() => new this.cls(n), 4198));
      return this.get(n, clear);
    }

    let ret = this.pools.get(n).next();
    if (clear) {
      let f = this.clearval;
      for (let i = 0; i < ret.length; i++) {
        ret[i] = f;
      }
    }

    return ret;
  }
}

let temps = new ArrayPool(Array, undefined);
let dtemps = new ArrayPool(Float64Array, 0.0);
let itemps = new ArrayPool(Int32Array, 0);

export function tempArray(n, clear) {
  return temps.get(n, clear);
}

export function dtempArray(n, clear) {
  return dtemps.get(n, clear);
}

export function itempArray(n, clear) {
  return itemps.get(n, clear);
}


export function svd_doubleArray(n, clear, name) {
  return dtempArray(n, clear);
}

export function svd_longArray(n, clear, name) {
  return itemparray(n, clear);
}

export function SAFE_FREE(f) {
}


export function svdWriteDenseArray() {

}

/**************************************************************
 * returns |a| if b is positive; else fsign returns -|a|      *
 **************************************************************/
export function svd_fsign(a, b) {
  if ((a >= 0.0 && b >= 0.0) || (a < 0.0 && b < 0.0))
    return (a);
  else
    return -a;
}

/**************************************************************
 * Function scales a vector by a constant.              *
 * Based on Fortran-77 routine from Linpack by J. Dongarra    *
 **************************************************************/
export function svd_dscal(n, da, dx, incx, start = 0) {
  let i;

  let di = start;

  if (n <= 0 || incx === 0) return;
  if (incx < 0) di += (-n + 1)*incx;

  for (i = 0; i < n; i++) {
    dx[di] *= da;
    di += incx;
  }
  return;
}

/**************************************************************
 * function scales a vector by a constant.              *
 * Based on Fortran-77 routine from Linpack by J. Dongarra    *
 **************************************************************/
export function svd_datx(n, da, dx, incx, dy, incy, startx = 0, starty = 0) {
  let i;
  let ix = startx, iy = starty;

  if (n <= 0 || incx === 0 || incy === 0 || da === 0.0) return;
  if (incx === 1 && incy === 1) {
    for (i = 0; i < n; i++) {
      dy[iy++] = da*dx[ix++];
      //*dy++ = da * (*dx++);
    }
  } else {
    if (incx < 0) ix += (-n + 1)*incx;
    if (incy < 0) iy += (-n + 1)*incy;

    for (i = 0; i < n; i++) {
      //*dy = da * (*dx);
      dy[iy] = da*dx[ix];
      ix += incx;
      iy += incy;
    }
  }
}

/**************************************************************
 * Function copies a vector x to a vector y              *
 * Based on Fortran-77 routine from Linpack by J. Dongarra    *
 **************************************************************/
export function svd_dcopy(n, dx, incx, dy, incy, startx = 0, starty = 0) {
  if (n <= 0 || incx === 0 || incy === 0) return;

  let i, ix = startx, iy = starty;

  if (incx === 1 && incy === 1) {
    for (i = 0; i < n; i++) {
      dy[iy++] = dx[ix++];
    }
  } else {
    if (incx < 0) ix += (-n + 1)*incx;
    if (incy < 0) iy += (-n + 1)*incy;
    for (i = 0; i < n; i++) {
      dy[iy] = dx[ix];

      ix += incx;
      iy += incy;
    }
  }
}

/**************************************************************
 * Function forms the dot product of two vectors.              *
 * Based on Fortran-77 routine from Linpack by J. Dongarra    *
 **************************************************************/
export function svd_ddot(n, dx, incx, dy, incy, startx = 0, starty = 0) {
  let i;
  let dot_product;

  if (n <= 0 || incx === 0 || incy === 0)
    return (0.0);

  let ix = startx, iy = starty;

  dot_product = 0.0;
  if (incx === 1 && incy === 1)
    for (i = 0; i < n; i++) {
      dot_product += dx[ix++]*dy[iy++];
    }
  else {
    if (incx < 0) ix += (-n + 1)*incx;
    if (incy < 0) iy += (-n + 1)*incy;
    for (i = 0; i < n; i++) {
      dot_product += dx[ix]*dy[iy]; //(*dx) * (*dy);
      ix += incx;
      iy += incy;
    }
  }

  return dot_product;
}

/**************************************************************
 * Constant times a vector plus a vector              *
 * Based on Fortran-77 routine from Linpack by J. Dongarra    *
 **************************************************************/
//export function svd_daxpy (long n, double da, double *dx, long incx, double *dy, long incy) {
export function svd_daxpy(n, da, dx, incx, dy, incy, startx = 0, starty = 0) {
  let i, ix = startx, iy = starty;

  if (n <= 0 || incx === 0 || incy === 0 || da === 0.0)
    return;

  if (incx === 1 && incy === 1) {
    for (i = 0; i < n; i++) {
      dy[iy] += da*dx[ix++];
      iy++;
    }
  } else {
    if (incx < 0) ix += (-n + 1)*incx;
    if (incy < 0) iy += (-n + 1)*incy;
    for (i = 0; i < n; i++) {
      dy[iy] += da*dx[ix];

      ix += incx;
      iy += incy;
    }
  }
  return;
}

/*********************************************************************
 * Function sorts array1 and array2 into increasing order for array1 *
 *********************************************************************/
//export function svd_dsort2(long igap, long n, double *array1, double *array2) {
export function svd_dsort2(igap, n, array1, array2) {
  let temp;
  let i, j, index;

  if (!igap)
    return;

  for (i = igap; i < n; i++) {
    j = i - igap;
    index = i;

    while (j >= 0 && array1[j] > array1[index]) {
      temp = array1[j];
      array1[j] = array1[index];
      array1[index] = temp;
      temp = array2[j];
      array2[j] = array2[index];
      array2[index] = temp;
      j -= igap;
      index = j + igap;
    }
  }

  svd_dsort2(igap/2, n, array1, array2);
}

/**************************************************************
 * Function interchanges two vectors                *
 * Based on Fortran-77 routine from Linpack by J. Dongarra    *
 **************************************************************/
//export function svd_dswap(long n, double *dx, long incx, double *dy, long incy) {
export function svd_dswap(n, dx, incx, dy, incy, startx = 0, starty = 0) {
  let i, ix = startx, iy = starty, dtemp;

  if (n <= 0 || incx === 0 || incy === 0)
    return;

  if (incx === 1 && incy === 1) {
    for (i = 0; i < n; i++) {
      dtemp = dy[iy];
      dy[iy++] = dx[ix];
      dx[ix++] = dtemp;
    }
  } else {
    if (incx < 0) ix += (-n + 1)*incx;
    if (incy < 0) iy += (-n + 1)*incy;
    for (i = 0; i < n; i++) {
      dtemp = dy[iy];
      dy[iy] = dx[ix];
      dx[ix] = dtemp;

      ix += incx;
      iy += incy;
    }
  }
}

/*****************************************************************
 * Function finds the index of element having max. absolute value*
 * based on FORTRAN 77 routine from Linpack by J. Dongarra       *
 *****************************************************************/
export function svd_idamax(n, dx, incx) {
  let ix, i, imax;
  let dtemp, dmax;

  let ix2;

  if (n < 1) return (-1);
  if (n === 1) return (0);
  if (incx === 0) return (-1);

  if (incx < 0) ix = (-n + 1)*incx;
  else ix = 0;

  imax = ix;
  ix2 += ix;
  dmax = fabs(dx[ix2]);

  for (i = 1; i < n; i++) {
    ix += incx;
    ix2 += incx;
    dtemp = fabs(dx[ix2]);
    if (dtemp > dmax) {
      dmax = dtemp;
      imax = ix;
    }
  }
  return (imax);
}

/**************************************************************
 * multiplication of matrix B by vector x, where B = A'A,     *
 * and A is nrow by ncol (nrow >> ncol). Hence, B is of order *
 * n = ncol (y stores product vector).                  *
 **************************************************************/
//export function svd_opb(SMat A, double *x, double *y, double *temp) {
export function svd_opb(A, x, y, temp) {
  let i, j, end;
  let pointr = A.pointr, rowind = A.rowind;
  let value = A.value;
  let n = A.cols;

  //SVDCount[SVD_MXV] += 2;
  for (let i = 0; i < n; i++) {
    y[i] = 0.0;
  }
  for (i = 0; i < A.rows; i++) {
    temp[i] = 0.0;
  }

  let yi = 0, xi = 0;

  for (i = 0; i < A.cols; i++) {
    end = pointr[i + 1];
    for (j = pointr[i]; j < end; j++) {
      temp[rowind[j]] += value[j]*x[xi];
    }
    xi++;
  }
  for (i = 0; i < A.cols; i++) {
    end = pointr[i + 1];

    for (j = pointr[i]; j < end; j++) {
      y[yi] += value[j]*temp[rowind[j]];
    }
    yi++;
  }
}

/***********************************************************
 * multiplication of matrix A by vector x, where A is     *
 * nrow by ncol (nrow >> ncol).  y stores product vector.  *
 ***********************************************************/
//export function svd_opa(SMat A, double *x, double *y) {
export function svd_opa(A, x, y) {
  let end, i, j;
  let pointr = A.pointr, rowind = A.rowind;
  let value = A.value;

  //SVDCount[SVD_MXV]++;
  for (let i = 0; i < A.rows; i++) {
    y[i] = 0.0;
  }

  for (i = 0; i < A.cols; i++) {
    end = pointr[i + 1];
    for (j = pointr[i]; j < end; j++) {
      y[rowind[j]] += value[j]*x[i];
    }
  }
}

/**************************************************************
 *                    *
 * Function finds sqrt(a^2 + b^2) without overflow or         *
 * destructive underflow.              *
 *                    *
 **************************************************************/
/**************************************************************

 Funtions used
 -------------

 UTILITY  dmax, dmin

 **************************************************************/
export function svd_pythag(a, b) {
  let p, r, s, t, u, temp;

  p = Math.max(fabs(a), fabs(b));
  if (p !== 0.0) {
    temp = Math.min(fabs(a), fabs(b))/p;
    r = temp*temp;
    t = 4.0 + r;
    while (t !== 4.0) {
      s = r/t;
      u = 1.0 + 2.0*s;
      p *= u;
      temp = s/u;
      r *= temp*temp;
      t = 4.0 + r;
    }
  }
  return (p);
}

