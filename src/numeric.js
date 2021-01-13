import * as numeric from './numeric.js';
import * as svd from './svd.js';

let exports = {};

for (let k in numeric) {
  exports[k] = numeric[k];
}

for (let k in numeric.generated) {
  exports[k] = numeric.generated[k];
}

for (let k in svd) {
  exports[k] = svd[k];
}

export default exports;
