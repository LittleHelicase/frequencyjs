
var _ = require("lodash");

var spectrumCreate = function(spectrumArray, sampling, signalLength){
  var spectrum = _.map(_.range(spectrumArray.length), function(idx){
    // frequency in hertz
    var freq = idx/signalLength * sampling;
    var freqFac = (idx == 0) ? 0.5 : 1;
    return { frequency: freq, amplitude: spectrumArray[idx]*freqFac };
  });
  spectrum.dominantFrequency = function() {
    return _.max(spectrum, function (s) {
      return s.amplitude;
    });
  };
  spectrum.amplitudes = function(){
    return _.map(spectrum, function(v){
      return v.amplitude;
    });
  }
  spectrum.sampling = sampling;
  return spectrum;
}

var Spectrum = function(spectrumArray, options){
  if("amplitudes" in spectrumArray)
    return spectrumArray;
  return spectrumCreate(spectrumArray, options.sampling, options.signalLength);
}

module.exports = Spectrum;
