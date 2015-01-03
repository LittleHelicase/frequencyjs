
module.exports = function(number){
  // bit check. every power of two has only one 1 in its binary representation
  // e.g. 32 = 100000
  // and a power of 2 minus 1 has all the lower bits set
  // e.g. 31 = 011111
  // the following statement uses a binary AND and it is always false for a 
  // non power of two, and true for a power of two
  // (it's probably the fastes way to do it)
  return ((number!=0) && !(number & (number-1)));
}
