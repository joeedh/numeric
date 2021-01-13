import numeric from './numeric.js';

let rand = new numeric.MersenneRandom();

let n = 4;
let mat = numeric.allocidentity(n, n);

for (let i=0; i<n; i++) {
  for (let j=0; j<n; j++) {
    mat[i][j] = rand.random();
  }
}

let v = new Array(n);
for (let i=0; i<n; i++) {
  v[i] = [rand.random()];
}

let eig=0;

for (let step=0; step<15; step++) {
  let v2 = numeric.dot(mat, v);
  console.log("");
  console.log(eig.toFixed(4));
  console.log(numeric.prettyPrint(v2))

  let tot = 0.0;

  for (let i = 0; i < n; i++) {
    tot += v2[i][0]**2;
    v[i][0] = v2[i][0];
  }

  tot = eig = Math.sqrt(tot);

  if (tot > 0.0) {
    tot = 1.0 / tot;

    for (let i=0; i<n; i++) {
      v[i][0] *= tot;
    }
  }
}

let ret = numeric.eig(mat);
console.log(numeric.prettyPrint(ret));
//console.log(numeric.neg)
//console.log(numeric.prettyPrint(mat));
