import numeric from './numeric.js';

window.numeric = numeric;

Math.fract = (f) => f - Math.floor(f);
Math.tent = (f) => 1.0 - Math.abs(Math.fract(f) - 0.5)*2.0;

function basicTests() {
  let mat = [
    [1, 0, 0.5, 0],
    [0, 1, 0.5, 0],
    [0.5, 0, 5, 0],
    [1, 2, 3, 1]
  ];

  console.log("DET", numeric.det(mat));
  let imat = numeric.pinv(mat);
  let mat2 = numeric.dot(mat, imat);

  console.log(numeric.print(mat2));
}
class Constraint {
  constructor(func, params, klst) {
    this.klst = [];
    this.glst = [];

    for (let i = 0; i < klst.length; i++) {
      this.klst.push(klst[i]);
      this.glst.push(new Float64Array(klst[i].length));
    }

    this.df = 0.00001;
    this.limit = 0.0001;

    this.func = func;
    this.params = params;
    this.k = 1.0;
  }

  evaluate(no_dvs = false) {
    let r1 = this.func(this.params);

    if (no_dvs) {
      return r1;
    }

    let klst = this.klst, glst = this.glst;
    let df = this.df;

    for (let i = 0; i < klst.length; i++) {
      let ks = klst[i], gs = glst[i];

      for (let j = 0; j < ks.length; j++) {
        let orig = ks[j];
        ks[j] += df;
        let r2 = this.func(this.params);
        ks[j] = orig;

        gs[j] = (r2 - r1)/df;
      }
    }

    return r1;
  }
}

class Solver {
  constructor() {
    this.cons = [];
    this.gk = 1.0;

    this.vec = undefined;
    this.klist = undefined;
    this.idxmap = undefined;
    this.ksmap = undefined;

    this.doinit = false;
  }

  add(con) {
    this.doinit = true;
    this.cons.push(con);
  }

  remove(con) {
    this.doinit = true;
    this.cons.remove(con);
  }

  solve_intern(gk, step, steps) {
    numeric.start();

    let vec = [];

    let klist = this.klist;
    let idxmap = this.idxmap;

    if (this.vec === undefined) {
      this.vec = new Float64Array(klist.length*2);
    }

    let cons = this.cons;

    let mat = new Array(cons.length);
    let err = 0.0;

    for (let i=0; i<mat.length; i++) {
      let row = mat[i] = new Float64Array(klist.length*2);

      for (let j=0; j<row.length; j++) {
        row[j] = 0.0;
      }
    }

    /*

    for (let con of cons) {
      let r1 = con.evaluate();

      if (Math.abs(r1) === 0.0) {
        continue;
      }

      let totg = 0.0;

      for (let i=0; i<con.klst.length; i++) {
        let ks = con.klst[i], gs = con.glst[i];

        for (let j=0; j<gs.length; j++) {
          totg += gs[j]*gs[j];
        }
      }

      if (totg === 0.0) {
        continue;
      }

      r1 /= totg;

      for (let i=0; i<con.klst.length; i++) {
        let ks = con.klst[i], gs = con.glst[i];
        let idx = this.idxmap.get(ks);

        for (let j=0; j<gs.length; j++) {
          let delta = -r1*gs[j]*con.k*gk;

          this.vec[idx+j] += delta;
          ks[j] += delta;
        }
      }
    }

    return;
    //*/

    for (let i = 0; i < cons.length; i++) {
      let con = this.cons[i];
      let row = mat[i];

      let r1 = con.evaluate();
      err += Math.abs(r1);

      vec.push(r1);

      for (let j=0; j<con.klst.length; j++) {
        let ks = con.klst[j], gs = con.glst[j];

        let idx = idxmap.get(ks);

        let mul = 20.0;

        row[idx] = gs[0]*mul;
        row[idx+1] = gs[1]*mul;
      }
    }

    for (let i=0; i<mat.length; i++) {
      let row = mat[i];
      for (let j=0; j<row.length; j++) {
        if (Math.abs(row[j]) > 0.000001) {
          row.start = j;
          break;
        }
      }
    }

    let perm = new Array(mat.length);
    let revperm = new Array(mat.length);
    for (let i=0; i<perm.length; i++) {
      perm[i] = i;
    }

    //perm.sort((a, b) => mat[a].start - mat[b].start);

    for (let i=0; i<perm.length; i++) {
      revperm[perm[i]] = i;
    }

    let mat2 = new Array(mat.length);
    let vec3 = new Array(vec.length);

    for (let i=0; i<mat2.length; i++) {
      mat2[i] = mat[perm[i]];
      vec3[i] = vec[perm[i]];
    }
    mat = mat2;

    if (step === 0) {
      _appstate.image = numeric.plotMatrix(mat);
    }

    let ret = numeric.leastSquares(mat, vec3);

    //console.log(ret);
    console.log(ret.length, klist.length*2);
    //return;

    let vec2 = this.vec;

    for (let i=0; i<vec3.length; i++) {
      let f = ret[perm[i]];

      if (isNaN(f) || !isFinite(f)) {
        continue;
      }

      let i2 = i>>1;
      let j = i & 1;

      let ks = klist[i2];
      ks[j] += -f;

      vec2[i] += -f;
    }

    numeric.end();

    console.log("error", err.toFixed(5));
    return err;
  }

  init() {
    this.doinit = false;

    let klist = this.klist = [];
    let idxmap = this.idxmap = new Map();
    let ksmap = this.ksmap = new Map();

    let visit = new WeakSet();
    for (let con of this.cons) {
      for (let ks of con.klst) {
        if (!visit.has(ks)) {
          idxmap.set(ks, klist.length*2);
          ksmap.set(klist.length*2, ks);
          klist.push(ks);
          visit.add(ks);
        }
      }
    }
  }

  solve(steps = 5, gk = 1.0) {
    gk *= this.gk;

    if (this.doinit) {
      this.init();
    }
    this.vec = undefined;

    for (let i = 0; i < steps; i++) {
      this.solve_intern(gk, i, steps);
    }

    return {
      vec : this.vec,
      ksmap : this.ksmap
    };
  }
}

function time_ms() {
  return performance.now();
}

class Vector2 {
  constructor(co, farray) {
    if (farray) {
      this.co = farray;
    } else {
      this.co = new Float64Array(2);
    }

    if (co instanceof Point) {
      co = co.co;
    }

    if (co) {
      this.co[0] = co[0];
      this.co[1] = co[1];
    }
  }

  get 0() {
    throw new Error("tried to access 0");
  }

  get 1() {
    throw new Error("tried to access 1");
  }

  get 2() {
    throw new Error("tried to access 2");
  }

  load(b) {
    if (b instanceof Vector2) {
      b = b.co;
    }

    this.co[0] = b[0];
    this.co[1] = b[1];

    return this;
  }

  zero() {
    this.co[0] = this.co[1] = 0.0;
    return this;
  }

  negate() {
    this.co[0] = -this.co[0];
    this.co[1] = -this.co[1];

    return this;
  }

  dot(b) {
    return this.co[0]*b.co[0] + this.co[1]*b.co[1];
  }

  vectorLength() {
    return Math.sqrt(this.dot(this));
  }

  vectorDistance(b) {
    let dx = this.co[0] - b.co[0];
    let dy = this.co[1] - b.co[1];

    return Math.sqrt(dx*dx + dy*dy);
  }

  normalize() {
    let len = this.co[0]*this.co[0] + this.co[1]*this.co[1];

    if (len < 0.00001) {
      return;
    }

    len = 1.0/Math.sqrt(len);

    this.co[0] *= len;
    this.co[1] *= len;

    return this;
  }

  add(b) {
    this.co[0] += b.co[0];
    this.co[1] += b.co[1];

    return this;
  }

  sub(b) {
    this.co[0] -= b.co[0];
    this.co[1] -= b.co[1];

    return this;
  }

  mul(b) {
    this.co[0] *= b.co[0];
    this.co[1] *= b.co[1];

    return this;
  }

  div(b) {
    this.co[0] /= b.co[0];
    this.co[1] /= b.co[1];

    return this;
  }

  mulScalar(b) {
    this.co[0] *= b;
    this.co[1] *= b;
    return this;
  }

  addScalar(b) {
    this.co[0] += b;
    this.co[1] += b;
    return this;
  }

  addFac(b, fac) {
    this.co[0] += b.co[0]*fac;
    this.co[1] += b.co[1]*fac;
    return this;
  }

  interp(b, t) {
    this.co[0] += (b.co[0] - this.co[0])*t;
    this.co[1] += (b.co[1] - this.co[1])*t;
    return this;
  }
}

class Point extends Vector2 {
  constructor(co, farray) {
    super(co, farray);

    this.startco = new Vector2(co);
    this.force = new Vector2();
    this.vel = new Vector2();
    this.oldco = new Vector2();
    this.w = 1.0;
    this.pin = 0.0;

    this.edges = [];
  }
}

class Edge {
  constructor(v1, v2) {
    this.v1 = v1;
    this.v2 = v2;

    this.length = 0;
  }

  otherVertex(v) {
    if (v === this.v1) {
      return this.v2;
    } else if (v === this.v2) {
      return this.v1;
    } else {
      throw new Error("vertex not in edge");
    }
  }
}

class AppState {
  constructor() {
    this.image = undefined;
    this.canvas = undefined;
    this.g = undefined;

    this.points = [];
    this.edges = [];
    this.comap = new Map();

    this.last_time = time_ms();
  }

  makeEdge(v1, v2) {
    for (let e of v1.edges) {
      if (e.otherVertex(v1) === v2) {
        return e;
      }
    }

    let e = new Edge(v1, v2);

    e.length = v1.vectorDistance(v2);

    v1.edges.push(e);
    v2.edges.push(e);

    this.edges.push(e);

    return e;
  }

  makePoint(co, fp) {
    let p = new Point(co, fp);
    this.points.push(p);
    this.comap.set(p.co, p);

    return p;
  }

  reset() {
    basicTests();

    this.comap = new Map();

    this.edges.length = 0;

    let ps = this.points;
    ps.length = 0;

    let totpoint = 32;
    let wid = Math.ceil(Math.sqrt(totpoint));

    totpoint = wid*wid;
    let parray = this.parray = new Float64Array(totpoint*2);
    let r = 0.35/wid;

    for (let i = 0; i < totpoint; i++) {
      let ix = i%wid, iy = ~~(i/wid);

      let x = (ix + 0.5)/wid, y = (iy + 0.5)/wid;

      x += (Math.random() - 0.5)*r;
      y += (Math.random() - 0.5)*r;

      let fp = new Float64Array(parray.buffer, i*2*8, 2);
      fp[0] = x;
      fp[1] = y;

      let p = this.makePoint(fp, fp);

      if (iy === wid-1) {
        p.pin = 1.0;
      }
    }

    function getp(x, y) {
      let idx = y*wid + x;

      return ps[idx];
    }

    for (let i = 0; i < wid - 1; i++) {
      for (let j = 0; j < wid - 1; j++) {
        let v1 = getp(i, j);
        let v2 = getp(i, j + 1);
        let v3 = getp(i + 1, j + 1);
        let v4 = getp(i + 1, j);

        this.makeEdge(v1, v2);
        this.makeEdge(v2, v3);
        this.makeEdge(v3, v4);
        this.makeEdge(v4, v1);
      }
    }
    return this;
  }

  start() {
    document.body.style.margin = "0px";
    document.body.style.padding = "0px";

    this.canvas = document.createElement("canvas");
    this.g = this.canvas.getContext("2d");

    document.body.appendChild(this.canvas);
    window.setInterval(this.on_tick.bind(this), 64);

    this.reset();
  }

  updateSize() {
    let canvas = this.canvas;
    let dpi = devicePixelRatio;

    let w = ~~((window.innerWidth - 25)*dpi);
    let h = ~~((window.innerHeight - 25)*dpi);

    if (canvas.width === w && canvas.height === h) {
      return false;
    }

    console.log("Size change detected");

    canvas.width = w;
    canvas.height = h;

    canvas.style["width"] = (w/dpi) + "px";
    canvas.style["height"] = (h/dpi) + "px";

    return true;
  }

  on_tick() {
    if (this.updateSize()) {
      this.draw();
    } else if (time_ms() - this.last_time > 45) {
      if (0) {
        this.step();
        window.redraw_all();
      }
      this.last_time = time_ms();
    }
  }

  step() {
    let gravity = 0.015;
    let damp = 0.999;

    for (let p of this.points) {
      p.force.zero();

      p.force.co[1] -= gravity;
      p.vel.add(p.force);
      //p.vel.mulScalar(damp);

      p.add(p.vel);
    }

    function clamp(f) {
      return Math.min(Math.max(f, 0.0), 1.0);
    }

    for (let p of this.points) {
      let x = clamp(p.co[0]);
      let y = clamp(p.co[1]);

      let dx = x - (p.co[0] + p.vel.co[0]);
      let dy = y - (p.co[1] + p.vel.co[1]);

      p.vel.co[0] += dx;
      p.vel.co[1] += dy;

      //p.co[0] = x;
      //p.co[1] = y;
    }

    this.solve();
  }

  solve() {
    let solver = new Solver();

    function len_c(params) {
      let v1 = params[0];
      let v2 = params[1];
      let restlen = params[2];

      return (v1.vectorDistance(v2) - restlen)**2;
    }

    for (let e of this.edges) {
      let con = new Constraint(len_c, [e.v1, e.v2, e.length], [e.v1.co, e.v2.co]);
      con.k = 0.1;
      solver.add(con);
    }

    for (let p of this.points) {
      if (p.pin === 0.0) {
        continue;
      }

      let con = new Constraint(len_c, [p, p.startco, 0.0], [p.co, p.startco.co]);
      solver.add(con);
    }

    for (let p of this.points) {
      p.oldco.load(p);
    }

    let {vec, ksmap} = solver.solve();

    for (let p of this.points) {
      p.load(p.oldco);
    }

    if (!vec) {
      return;
    }

    let comap = this.comap;
    for (let i=0; i<vec.length; i += 2) {
      let p = comap.get(ksmap.get(i));

      if (!p) {
        continue;
      }
      p.vel.co[0] += vec[i];
      p.vel.co[1] += vec[i+1];
    }
  }

  on_keydown(e) {
    console.log(e.keyCode);
    switch (e.keyCode) {
      case 68: //dkey
        this.step();
        window.redraw_all();
        break;
      case 82: //rkey
        this.reset();
        window.redraw_all();
        break;
    }
  }

  draw() {
    this.step();

    console.log("Draw");
    this.updateSize();

    let canvas = this.canvas, g = this.g;

    g.clearRect(0, 0, canvas.width, canvas.height);

    if (this.image) {
      g.putImageData(this.image, 0, 0);
    }

    g.save();

    g.strokeStyle = "black";

    let scale1 = Math.min(canvas.width, canvas.height);
    let scale2 = 0.75;
    let scale = scale1*scale2;

    g.lineWidth /= scale1;

    //scale and flip y
    g.scale(scale1, -scale1);
    g.translate(0.0, -scale2);
    g.scale(scale2, scale2);
    g.translate(0.1, -0.1);

    g.beginPath();
    g.rect(0, 0, 1.0, 1.0);
    g.stroke();

    g.beginPath();

    for (let e of this.edges) {
      g.moveTo(e.v1.co[0], e.v1.co[1]);
      g.lineTo(e.v2.co[0], e.v2.co[1]);
    }
    g.stroke();

    g.beginPath();
    let r = 0.1/Math.sqrt(this.points.length);

    for (let p of this.points) {
      g.moveTo(p.co[0], p.co[1]);
      g.arc(p.co[0], p.co[1], r, -Math.PI, Math.PI);
    }
    g.fill();

    g.restore();
    //window.redraw_all()
  }
}

export function start() {
  let animreq;

  function draw() {
    animreq = undefined;
    _appstate.draw();
  }

  window.redraw_all = function () {
    if (animreq) {
      return;
    }

    animreq = window.requestAnimationFrame(draw);
  }

  window.addEventListener("keydown", (e) => _appstate.on_keydown(e));

  window._appstate = new AppState();

  _appstate.start();
  window.redraw_all();
}
