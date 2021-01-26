export function plotMatrix(matrix) {
  let m = matrix.length;
  let n = matrix[0].length;

  let msize = Math.max(n, m);
  let cellsize;

  if (msize > 128) {
    cellsize = 2;
  } else if (msize > 64) {
    cellsize = 4;
  } else if (msize > 32) {
    cellsize = 8;
  }

  let width = cellsize*n;
  let height = cellsize*m;

  let image = new ImageData(width, height);
  let idata = image.data;

  for (let i=0; i<idata.length; i++) {
    idata[i] = 255;
  }

  let min=1e17, max=-1e17;
  for (let i=0; i<m; i++) {
    for (let j = 0; j < n; j++) {
      min = Math.min(min, matrix[i][j]);
      max = Math.max(max, matrix[i][j]);
    }
  }

  if (min === max) {
    return image;
  }

  for (let i=0; i<m; i++) {
    for (let j=0; j<n; j++) {
      let x = j*cellsize;
      let y = i*cellsize;

      let idx = (y*width + x)*4;
      let f = matrix[i][j];

      let f2 = (f - min) / (max - min);

      for (let y1=0; y1<cellsize; y1++) {
        for (let x1 = 0; x1 < cellsize; x1++) {
          let x2 = x + x1;
          let y2 = y + y1;
          if (x2 >= width || y2 >= height) {
            continue;
          }

          let idx = (y2*width + x2)*4;

          if (f > 0.0) {
            idata[idx] = ~~(f2*255);
          } else if (f < 0.0) {
            idata[idx + 1] = ~~(f2*255);
          }

          if (Math.abs(f) > 0.000001) {
            idata[idx + 2] = 255;
          }
        }
      }
    }
  }

  return image;
}