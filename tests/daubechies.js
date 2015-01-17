
var chai = require("chai");
var _ = require("lodash");
chai.should();

var Generator = require("../lib/generator.js");
var Transform = require("../lib/transform.js");
var Processing = require("../lib/processing.js");
var Signals = require("../lib/utils/signals.js");

describe("Daubechies Wavelet Transform", function(){
  it("should transform and backtransform an impuls correctly", function(){
    var impuls = [1,2,3,4,5,6,7,8];
    var spec = Transform
      .toSpectrum(impuls,{method:"daubechies"});
    var impBack = Transform.toSignal(spec,{method:"daubechies"});
    var areEqual = Processing.equal(impBack,impuls);
    areEqual.should.be.true;
  });
  it("should be invertible", function(){
    var signal = Generator.sines([
      {frequency: 40, amplitude: 0.5},
      {frequency: 12, amplitude: 0.8},
      ]).create({length:128, sampling: 100});
    var spec = Transform.toSpectrum(signal,{method:"daubechies"});
    var signalBack = Transform.toSignal(spec,{method:"daubechies"});
    var areEqual = Processing.equal(signal, signalBack);
    areEqual.should.be.true;
  });
});
