import * as numeric from './base.js';
import * as svd from './svd.js';
import * as plot from './plot.js';

let exports = {};

let direct = new Set([
  "precision", "epsilon", "_startdepth", "largeArray"
]);

for (let k of direct) {
  Object.defineProperty(exports, k, {
    get() {
      return numeric[k]
    },
    set(v) {
      numeric[k] = v;
    }
  });
}

for (let k in numeric) {
  if (direct.has(k)) {
    continue;
  }

  exports[k] = numeric[k];
}

for (let k in numeric.generated) {
  if (direct.has(k)) {
    continue;
  }

  exports[k] = numeric.generated[k];
}

for (let k in svd) {
  exports[k] = svd[k];
}

for (let k in plot) {
  exports[k] = plot[k];
}

export default exports;
