
module.exports = {
  triangle: function(period){
    return function(idx){
      return 1 - (Math.abs(((idx + period*0.25) % (period)) - period*0.5)) / period*4;
    }
  },
  rectangle: function(period){
    return function(idx){
      return (idx % period) < (period*0.5) ? 1 : -1;
    }
  }
}
