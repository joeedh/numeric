export const freetag = Symbol("freetag");
export const nexttag = Symbol("nexttag");
export const prevtag = Symbol("prevtag");
export const pooltag = Symbol("pooltag");
export const scopetag = Symbol("scopetag");

export var scope;

export class cachepool extends Array {
  constructor(cb, initialCount = 8) {
    super();

    this.factory = cb;

    this.freelist = [];

    for (let i = 0; i < initialCount; i++) {
      let item = this._alloc();
      item[freetag] = true;
      this.freelist.push(item);
    }
  }

  _alloc() {
    let ret = this.factory();

    ret[freetag] = false;
    ret[nexttag] = undefined;
    ret[prevtag] = undefined;
    ret[pooltag] = this;

    this.push(ret);

    return ret;
  }

  alloc() {
    if (this.freelist.length > 0) {
      let item = this.freelist.pop();

      item[freetag] = false;
      scope.addBlock(item);

      return item;
    }

    return this._alloc();
  }

  free(item) {
    if (item[freetag]) {
      console.warn("double free");
      return;
    }

    if (item[scopetag]) {
      item[scopetag].remBlock(item);
    }

    item[freetag] = true;
    this.freelist.push(item);

    return this;
  }

  freeAll() {
    for (let item of this) {
      if (!item[freetag]) {
        if (item[scopetag]) {
          item[scopetag].remBlock(item);
        }

        item[freetag] = true;
        this.freelist.push(item);
      }
    }

    return this;
  }
}

let arraypools = new Map();
let matpools = new Map();
let _startdepth = 0;

export const pools = [arraypools, matpools];

export function allocmat(m, n) {
  let key = m*8192 + n;
  //key = key | (1 < 23);

  //key = "m" + m + ":" + n;
  let pool = matpools.get(key);
  if (!pool) {
    pool = new cachepool(() => {
      let ret = new Array(m);

      for (let i = 0; i < m; i++) {
        ret.push(new Array(n));
      }

      return ret;
    }, 8);

    matpools.set(key, pool);
  }

  return pool.alloc();
}

export function allocarray(n) {
  let key = n;

  //key = "a:" + n;
  let pool = arraypools.get(key);

  if (!pool) {
    pool = new cachepool(() => {
      return new Array(n);
    }, 8);

    arraypools.set(key, pool);
  }

  return pool.alloc();
}

export function free(block) {
  block[pooltag].free(block);
}

export class Scope {
  constructor() {
    this.blocklist = {
      first: undefined,
      last : undefined
    };

    this.parent = undefined;
  }

  //adds block b and any subblocks it contains to this scope
  addBlock(b) {
    if (b[scopetag] === this) {
      return; //block is already in this scope
    }

    if (b[scopetag] !== undefined) {
      b[scopetag].remBlock(b);
    }

    b[scopetag] = this;

    if (Array.isArray(b)) {
      for (let block2 of b) {
        if (typeof block2 === "object" && block2[pooltag]) {
          this.addBlock(block2);
        }
      }
    }

    let bl = this.blocklist;

    if (!bl.first) {
      b[nexttag] = b[prevtag] = undefined;
      bl.first = bl.last = b;
      return;
    }

    b[nexttag] = undefined;
    b[prevtag] = bl.last;
    bl.last[nexttag] = b;
    bl.last = b;
  }

  reset() {
    this.blocklist.first = this.blocklist.last = undefined;
    this.parent = undefined;

    return this;
  }

  remBlock(b) {
    if (b[scopetag] !== this) {
      throw new Error("block not part of scope");
    }

    let bl = this.blocklist;

    if (b === bl.first) {
      bl.first = b[nexttag];
    }
    if (b === gl.last) {
      bl.last = b[prevtag];
    }

    if (b[prevtag]) {
      b[prevtag][nexttag] = b[nexttag];
    }

    if (b[nexttag]) {
      b[nexttag][prevtav] = b[prevtab];
    }

    b[scopetag] = undefined;

    return b;
  }

  end() {
    let next;

    for (let b = this.blocklist.first; b; b = next) {
      next = b[nexttag];

      if (!b[freetag]) {
        b[pooltag].free(b);
      }
    }
  }
}

scope = new Scope();

let scopealloc = new cachepool(() => new Scope(), 8);

export function start() {
  if (_startdepth === 0) {
    for (let pools of pools) {
      for (let pool of pools.values()) {
        pool.freeAll();
      }
    }
  }

  let scope2 = scopealloc.alloc().reset();

  scope2.parent = scope;
  scope = scope2;

  _startdepth++;
}

export function scopePull(block) {
  if (scope.parent && block[scopetag] !== scope.parent) {
    scope.parent.addBlock(block);
  }
}

export function bench() {
  return performance.now();
}

export function list(iter) {
  let ret = [];

  for (let item of iter) {
    ret.push(item);
  }

  return ret;
}

//arguments are list of blocks to be moved to parent scope
export function end() {
  for (let i = 0; i < arguments.length; i++) {
    let b = arguments[i];
    let pscope = scope.parent;

    if (pscope) {
      pscope.addBlock(b);
    }
  }

  _startdepth--;

  if (_startdepth < 0) {
    console.error("Mismatched numeric.end()");
    _startdepth = 0;
    return;
  }

  let oldscope = scope;

  scope = scope.parent;

  oldscope.end();
  scopealloc.free(oldscope);
}

// 1. Utility function

export class cachering extends Array {
  constructor(cb, count) {
    super();

    this.length = count;
    this.cur = 0;

    for (let i = 0; i < this.length; i++) {
      this[i] = cb();
    }
  }

  next() {
    let ret = this[this.cur];
    this.cur = (this.cur + 1)%this.length;

    return ret;
  }
}

export function allocidentity(n) {
  let mat = allocmat(n, n);

  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      mat[i][j] = 0.0;
    }
  }

  for (let i = 0; i < n; i++) {
    mat[i][i] = 1.0;
  }

  return mat;
}

export function allocclone(A) {
  let m = A.length, n = A[0].length;

  let mat = allocmat(m, n);

  for (let i = 0; i < m; i++) {
    let row1 = mat[i];
    let row2 = A[i];

    for (let j = 0; j < n; j++) {
      row1[j] = row2[j];
    }
  }

  return mat;
}

export function _myIndexOf (w) {
  var n = this.length,k;
  for(k=0;k<n;++k) if(this[k]===w) return k;
  return -1;
};
export const myIndexOf = (Array.prototype.indexOf)?Array.prototype.indexOf:numeric._myIndexOf;
