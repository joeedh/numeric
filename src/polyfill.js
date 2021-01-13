if (Array.prototype.remove === undefined) {
  Array.prototype.remove = function (item, suppress_error) {
    var i = this.indexOf(item);

    if (i < 0) {
      if (suppress_error)
        console.trace("Warning: item not in array", item);
      else
        throw new Error("Error: item not in array " + item);

      return;
    }

    while (i < this.length-1) {
      this[i] = this[i+1];
      i++;
    }

    this.length--;

    return this;
  }
}
