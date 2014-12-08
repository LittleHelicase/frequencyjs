
var chai = require("chai");
chai.should();

var Signal = require("../lib/signal.js");
var _ = require("lodash");

describe("Signal Equality", function(){
  it("a signal should be equal to itself", function(){
    var sig = _.range(27);
    var areEqual = Signal.equal(sig, sig);
    areEqual.should.be.true;
  });
  it("a signal should not equal a shorter version of itself", function(){
    var sig = _.range(29);
    var sig2 = sig.slice(0,27);
    var areEqual = Signal.equal(sig,sig2);
    areEqual.should.be.false;
  });
  it("a signal should not equal a similar version of itself", function(){
    var sig = _.range(17);
    var sig2 = _.map(sig, function(v){ return v + 1E-2});
    var areEqual = Signal.equal(sig,sig2);
    areEqual.should.be.false;
  });
  it("should not fail for long signals", function(){
    var sig = _.range(10);
    var sig2 = sig.slice();
    sig2[0]=0.1;
    var areEqual = Signal.equal(sig,sig2);
    areEqual.should.be.false;
  })
});
