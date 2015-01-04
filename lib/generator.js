/**
  * Generators for different signals (currently only sine waves)
  */

// for now we use dsp.js
var dsp = require("../dsp.js");
var _ = require("lodash");
var Signal = require("./signal.js");

/** A generator can generate a periodic function
 *  using a given function
 */
var generator = function(genFunction, cfg){
  
  var data = []
  return {
    /** creates a signal
     *  call
     *   - create() [uses default: sampling=44100, length=1000]
     *   - create({sampling: 44100, length: 1000})
     */
    create: function(options){
      options = _.defaults(options || {}, {sampling:44100, length:1000});
//      var generator = new dsp.Oscillator(cfg.__dspGenerator, 
//        cfg.frequency, cfg.amplitude, options.length, options.sampling);
//      generator.generate();
      var signalData = _.map(_.range(options.length), function(idx){
        return cfg.amplitude * genFunction(cfg.functionPeriod * cfg.frequency * idx / options.sampling);
      });
//      var signalData = generator.signal;
      var s = Signal(signalData, options);
      s.amplitude = cfg.amplitude;
      s.frequency = cfg.frequency;
      return s;
    }
  }
}

var mergeTwoSignals = function(s1,s2){
  var sigData = _(s1)
    .zip(s2)
    .map(function(p){
      return p[0].value + p[1].value;
    })
    .value();
  return Signal(sigData,{sampling:s1.sampling});
}

var mergeSignals = function(sigs){
  if(_.isArray(sigs)){
    return _.reduce(sigs,mergeTwoSignals);
  } else {
    return _.reduce(arguments,mergeTwoSignals);
  }
}

var combinator = function(generators){
  return {
    create: function(options){
      options = _.defaults(options || {}, {sampling:44100, length:1000});
      return _(generators)
        .invoke("create",options)
        .reduce(mergeTwoSignals);
    }
  }
}

var Signals = require("./utils/signals.js");

var Generators = {
  /** sine wave generator.
   *  call
   *   - sine() [uses default: frequency=440, phase = 0, amplitude: 1]
   *   - sine({frequency:440, phase:0, amplitude: 1})
   *  returns a generator for the sine wave
   */
  sine: function(options){
    options = _.defaults(options || {}, {frequency: 440, phase: 0, amplitude: 1});
    options.functionPeriod = 2*Math.PI;
    return generator(Math.sin, options);
  },
  sines: function(signals){
    var gens = _.map(signals, function(s){
      return Generators.sine(s); });
    return combinator(gens);
  },
  triangle: function(options){
    options = _.defaults(options || {}, {frequency: 440, phase: 0, amplitude: 1});
    options.functionPeriod = 1;
    return generator(Signals.triangle(1), options);
  },
  rectangle: function(options){
    options = _.defaults(options || {}, {frequency: 440, phase: 0, amplitude: 1});
    options.functionPeriod = 1;
    return generator(Signals.rectangle(1), options);
  }
}

module.exports = Generators;
