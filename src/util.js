export const freetag = Symbol("freetag");
export const nexttag = Symbol("nexttag");
export const prevtag = Symbol("prevtag");
export const pooltag = Symbol("pooltag");
export const scopetag = Symbol("scopetag");

export var scope;

let arraypools = new Map();
let matpools = new Map();
export var _startdepth = 0;

export const pools = [arraypools, matpools];


//from:https://en.wikipedia.org/wiki/Mersenne_Twister
function _int32(x) {
  // Get the 31 least significant bits.
  return ~~(((1<<30) - 1) & (~~x))
}

export class MersenneRandom {
  constructor(seed) {
    // Initialize the index to 0
    this.index = 624;
    this.mt = new Uint32Array(624);

    this.seed(seed);
  }

  random() {
    return this.extract_number()/(1<<30);
  }

  seed(seed) {
    seed = ~~(seed*8192);

    // Initialize the index to 0
    this.index = 624;
    this.mt.fill(0, 0, this.mt.length);

    this.mt[0] = seed;  // Initialize the initial state to the seed

    for (var i = 1; i < 624; i++) {
      this.mt[i] = _int32(
        1812433253*(this.mt[i - 1] ^ this.mt[i - 1]>>30) + i);
    }
  }

  extract_number() {
    if (this.index >= 624)
      this.twist();

    var y = this.mt[this.index];

    // Right shift by 11 bits
    y = y ^ y>>11;
    // Shift y left by 7 and take the bitwise and of 2636928640
    y = y ^ y<<7 & 2636928640;
    // Shift y left by 15 and take the bitwise and of y and 4022730752
    y = y ^ y<<15 & 4022730752;
    // Right shift by 18 bits
    y = y ^ y>>18;

    this.index = this.index + 1;

    return _int32(y);
  }

  twist() {
    for (var i = 0; i < 624; i++) {
      // Get the most significant bit and add it to the less significant
      // bits of the next number
      var y = _int32((this.mt[i] & 0x80000000) +
        (this.mt[(i + 1)%624] & 0x7fffffff));
      this.mt[i] = this.mt[(i + 397)%624] ^ y>>1;

      if (y%2 != 0)
        this.mt[i] = this.mt[i] ^ 0x9908b0df;
    }

    this.index = 0;
  }
}

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
    let item;

    if (this.freelist.length > 0) {
      item = this.freelist.pop();
    } else {
      item = this._alloc();
    }

    item[pooltag] = this;
    item[freetag] = false;
    scope.addBlock(item);

    return item;
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

export class Matrix extends Array {
  constructor() {
    super();
  }

  set 0(v) {
    if (typeof v === "string") {
      throw new Error("eek");
    }

    this._zero = v;
  }

  get 0() {
    return this._zero;
  }
}

export function allocmat(m, n=m) {
  let key = m*8192 + n;
  //key = key | (1 < 23);

  //key = "m" + m + ":" + n;
  let pool = matpools.get(key);
  if (!pool) {
    pool = new cachepool(() => {
      let ret = new Matrix(m);

      for (let i = 0; i < m; i++) {
        ret[i] = new Array(n);
      }

      return ret;
    }, 8);

    matpools.set(key, pool);
  }

  let ret = pool.alloc();
  ret.length = m;

  for (let i=0; i<ret.length; i++) {
    if (typeof ret[i] === "string") {
      console.warn("Corrupted matrix; fixing.", i);
      ret[i] = new Array(n);
    }

    ret[i].length = n;
  }

  return ret;
}

export function allocarray(n) {
  let key = n;

  if (isNaN(n)) {
    console.log(n);
    throw new Error("invalid array length" + n);
  }

  //key = "a:" + n;
  let pool = arraypools.get(key);

  if (!pool) {
    pool = new cachepool(() => {
      return new Array(n);
    }, 8);

    arraypools.set(key, pool);
  }

  let ret = pool.alloc();
  ret.length = n;

  return ret;
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
    if (b === bl.last) {
      bl.last = b[prevtag];
    }

    if (b[prevtag]) {
      b[prevtag][nexttag] = b[nexttag];
    }

    if (b[nexttag]) {
      b[nexttag][prevtag] = b[prevtag];
    }

    b[scopetag] = undefined;

    return b;
  }

  end() {
    let next;

    for (let b = this.blocklist.first; b; b = next) {
      next = b[nexttag];

      if (!b[freetag]) {
        if (!b[pooltag]) {
          console.warn("Missing pool tag", b.constructor);
          continue;
        }

        b[pooltag].free(b);
      }
    }
  }
}

scope = new Scope();

let scopealloc = new cachepool(() => new Scope(), 8);

export function _freeAll() {
  for (let poolmap of pools) {
    for (let pool of poolmap.values()) {
      pool.freeAll();
    }
  }
}

export function start() {
  if (_startdepth === 0) {
    _freeAll();
  }

  let scope2 = scopealloc.alloc().reset();

  scope2.parent = scope;
  scope = scope2;

  _startdepth++;
}

export function reportMemory() {
  let blocks = [];

  for (let poolset of pools) {
    for (let pool of poolset.values()) {
      for (let block of pool) {
        if (!block[freetag]) {
          blocks.push(block);
        }
      }
    }
  }

  let tot = 0;
  for (let block of blocks) {
    tot += block.length*16;

    console.log(block.constructor.name);
  }

  tot /= 1024;
  tot = tot.toFixed(2);
  console.log("estimated leaked memory:", tot + "kb");
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
    //console.error("Mismatched numeric.end()");
    throw new Error("Mismatched numeric.end()");

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
  let m = A.length;

  if (typeof A[0] === "number") {
    let ret = allocarray(A.length);

    for (let i=0; i<ret.length; i++) {
      ret[i] = A[i];
    }

    return ret;
  }

  let n = A[0].length;
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
