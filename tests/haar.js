
var chai = require("chai");
var _ = require("lodash");
chai.should();

var Generator = require("../lib/generator.js");
var Transform = require("../lib/transform.js");
var Processing = require("../lib/processing.js");
var Signals = require("../lib/utils/signals.js");

describe("Haar Wavelet Transform", function(){
  /*it("should detect the frequency of the sine generator", function(){
    // generate signal of 440 hertz
    var signal = Generator
    .sine({frequency: 440})
    .create({length: 100, sampling: 4400});
    
    // calculate the dominant frequency of the signal using dft
    var domFreq = Transform
    .toSpectrum(signal,{method: "haar"})
    .dominantFrequency();
    domFreq.frequency.should.equal(440);
  });*/
/*  it("should detect different frequencies", function(){
    //for(var i=1; i<4; i++) {
    for(var i=1; i<4; i++) {
      var rect = Signals.rectangle(Math.pow(2,4-i));
      cosFunc = _.map(_.range(8), function (idx) {
        return rect(idx);
      });
      console.log("func");
      console.log(cosFunc);
      // calculate the dominant frequency of the signal using dft
      var spec = Transform.toSpectrum(cosFunc, {method: "haar", sampling:8});
      console.log(spec);
      var domFreq = Transform
      .toSpectrum(cosFunc,{method: "haar", sampling:8})
      .dominantFrequency();
      domFreq.frequency.should.equal(i)
    }
  });*/
});
