
var chai = require("chai");
var _ = require("lodash");
chai.should();

var Generator = require("../lib/generator.js");
var Transform = require("../lib/transform.js");
var Processing = require("../lib/processing.js");
var Spectrum = require("../lib/spectrum.js");

describe("Haar Wavelet Transform", function(){
  it("should detect different (rectangle) frequencies", function(){
    //for(var i=1; i<4; i++) {
    for(var i=1; i<4; i++) {
      var freq = Math.pow(2,i-1);
      var rectFunc = Generator.rectangle({frequency:freq})
        .create({sampling:8,length:8});
      var spec = Transform
        .toSpectrum(rectFunc,{method: "haar"})
      // calculate the dominant frequency of the signal using dft
      var domFreq = Transform
      .toSpectrum(rectFunc,{method: "haar"})
      .dominantFrequency();
      domFreq.frequency.should.equal(Math.pow(2,i - 1))
    }
  });
  it("should transform an impuls correctly", function(){
    // example taken from http://en.wikipedia.org/wiki/Discrete_wavelet_transform#Comparison_with_Fourier_transform
    var impuls = [1,0,0,0];
    var spec = Transform
      .toSpectrum(impuls,{method:"haar", sampling: 4})
      .amplitudes();
    var areEqual = Processing.equal(spec, [0.25,0.25,0.5,0]);
    areEqual.should.be.true;
  });
  it("should be invertible", function(){
    var spec = Spectrum([0.25,0.25,0.5,0],{timeDependent: true, sampling:4});
    var signal = Transform
      .toSignal(spec,{method:"haar"});
    signal.sampling.should.be.equal(4);
    console.log(signal);
    var areEqual = Processing.equal(signal, [1,0,0,0]);
    areEqual.should.be.true;
/*    var signal = Generator.sines([
      {frequency: 40, amplitude: 0.5},
      {frequency: 12, amplitude: 0.8},
      ]).create({length:128, sampling: 100});
    var spec = Transform.toSpectrum(signal,{method:"haar"});
    var signalBack = Transform.toSignal(spec,{method:"haar"});
    var areEqual = Processing.equal(signal, signalBack);
    areEqual.should.be.true;*/
  });
});
