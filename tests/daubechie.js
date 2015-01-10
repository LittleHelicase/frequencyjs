
var chai = require("chai");
var _ = require("lodash");
chai.should();

var Generator = require("../lib/generator.js");
var Transform = require("../lib/transform.js");
var Processing = require("../lib/processing.js");
var Signals = require("../lib/utils/signals.js");

describe("Daubechies Wavelet Transform", function(){
  it("should equal haar for D2", function(){
    var sineSig = Generator.sine().create({length:128});
    var haarSpec = Transform.toSpectrum(sineSig, {method: "haar",sampling:8});
    var daubSpec = Transform.toSpectrum(sineSig, {method: "daubechies", taps:2,sampling:8});
    haarSpec.dominantFrequency().frequency.should.equal(daubSpec.dominantFrequency().frequency);
  });
});
