/**
  * Generators for different signals (currently only sine waves)
  */

// for now we use dsp.js
var dsp = require("../dsp.js");
var _ = require("lodash");

/** A generator can generate a periodic function
 *  using a given function
 */
var Generator = function(genFunction, cfg){
  
  var data = []
  return {
    /** creates a signal
     *  call
     *   - create() [uses default: sampling=44100, length=1000]
     *   - create({sampling: 44100, length: 1000})
     */
    create: function(options){
      options = _.defaults(options || {}, {sampling:44100, length:1000});
      var generator = new dsp.Oscillator(cfg.__dspGenerator, 
        cfg.frequency, cfg.amplitude, options.length, options.sampling);
      generator.generate();
      var signal = generator.signal;
      var curT = 0;
      data = _.map(signal, function(v){
        ret = {t: curT, value: v};
        curT = curT + 1/options.sampling;
        return ret;
      });
      data.amplitude = cfg.amplitude;
      data.frequency = cfg.frequency;
      data.sampling = options.sampling;
      data.dt = 1/options.sampling;
      return this;
    },
    /**
      * returns the values of the signal without the corresponding time. The
      * the index of the value indicates the time. At index i the value 
      * corresponds to time:
      *   t = i * dt = i / samplingRate
      */
    values: function(){
      return _.map(data, function(v){ return v.value; });
    },
    /**
      * returns an array of objects that contain the time and value at
      * the given time. The values are always ordered ascending with time
      *   [{t:0, value: 0}, {t:0.1, value: 0.023}, ... ]
      */
    data: function(){
      return data;
    }
  }
}

module.exports = {
  /** sine wave generator.
   *  call
   *   - sine() [uses default: frequency=440, phase = 0, amplitude: 1]
   *   - sine({frequency:440, phase:0, amplitude: 1})
   *  returns a generator for the sine wave
   */
  sine: function(options){
    options = _.defaults(options || {}, {frequency: 440, phase: 0, amplitude: 1});
    options.functionPeriod = 2*Math.PI;
    // only for use with dsp... gets replaces..
    options.__dspGenerator = dsp.SINE;
    return Generator(Math.sin, options);
  }
}