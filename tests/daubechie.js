
var chai = require("chai");
var _ = require("lodash");
chai.should();

var Generator = require("../lib/generator.js");
var Transform = require("../lib/transform.js");
var Processing = require("../lib/processing.js");
var Signals = require("../lib/utils/signals.js");

describe("Daubechies Wavelet Transform", function(){
  it("should equal haar for D2", function(){
    var sineSig = Generator.sine().create({length:1024});
    var haarSpec = Transform.toSpectrum(sineSig, {method: "haar",sampling:8});
    var daubSpec = Transform.toSpectrum(sineSig, {method: "daubechies", taps:2,sampling:8});
    haarSpec.dominantFrequency().frequency.should.equal(daubSpec.dominantFrequency().frequency);
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
