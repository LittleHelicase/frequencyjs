var charts = {
  "transform_signal" : {
    width: 500,
    height: 300,
    data: _.map(_.range(200), function(idx){ return {t: idx/ 200, value: Math.sin(4*idx*2*Math.PI/199) }; }),
    ticks: [0,0.5,1],
    type: "signal"
  },
  "transform_spectrum" : {
    width: 500,
    height: 300,
    data: _.map(_.range(10), function(idx){ return {t: idx, value: ((idx == 4) ? 1 : 0) }; }),
    ticks: _.range(10),
    type: "spectrum"
  }
};

// taken from https://gist.github.com/benjchristensen/2579599. Thx to benjchristensen
function basePlot(chart, key){
  var m = [40, 80, 40, 80]; // margins
  var w = chart.width - m[1] - m[3]; // width
  var h = chart.height - m[0] - m[2]; // height

  var data = chart.data;

  var x = d3.scale.linear().domain([0, d3.max(data, function(d){return d.t; })]).range([0, w]);
  var y = d3.scale.linear().domain([-1, 1]).range([h, 0]);
  // automatically determining max range can work something like this
  // var y = d3.scale.linear().domain([0, d3.max(data)]).range([h, 0]);

  var svg = d3.select("#"+key).append("svg:svg")
    .attr("width", w + m[1] + m[3])
    .attr("height", h + m[0] + m[2]);
  
  svg.append("rect")
    .attr("width", chart.width)
    .attr("height", chart.height);
  
  var graph = svg
    .append("svg:g")
    .attr("transform", "translate(" + m[3] + "," + m[0] + ")");

  var xAxis = d3.svg.axis().scale(x).tickValues(chart.ticks);
  graph.append("svg:g")
    .attr("class", "x axis")
    .attr("transform", "translate(0," + y(0) + ")")
    .call(xAxis);


  var yAxisLeft = d3.svg.axis().scale(y).ticks(4).orient("left");
  graph.append("svg:g")
    .attr("class", "y axis")
    .attr("transform", "translate(-25,0)")
    .call(yAxisLeft);
  return graph;
}

function plot(chart,key, graph){
  var m = [40, 80, 40, 80]; // margins
  var w = chart.width - m[1] - m[3]; // width
  var h = chart.height - m[0] - m[2]; // height

  var data = chart.data;

  var x = d3.scale.linear().domain([0, d3.max(data, function(d){return d.t; })]).range([0, w]);
  var y = d3.scale.linear().domain([-1, 1]).range([h, 0]);
  // automatically determining max range can work something like this
  // var y = d3.scale.linear().domain([0, d3.max(data)]).range([h, 0]);

  var line = d3.svg.line()
    .x(function(d,i) { 
        return x(d.t); 
    })
    .y(function(d) { 
        return y(d.value); 
    });

  graph.append("svg:path").attr("class","function").attr("d", line(data));
}

function plotSpectrum(chart,key,graph){
  var m = [40, 80, 40, 80]; // margins
  var w = chart.width - m[1] - m[3]; // width
  var h = chart.height - m[0] - m[2]; // height

  var data = _(chart.data).map(function(v){ 
    return [
      {t:v.t, value: 0 },
      {t:v.t, value: v.value },
      {t:v.t, value: 0 }]
    })
    .flatten()
    .value();

  var x = d3.scale.linear().domain([0, d3.max(data, function(d){return d.t; })]).range([0, w]);
  var y = d3.scale.linear().domain([-1, 1]).range([h, 0]);
  // automatically determining max range can work something like this
  // var y = d3.scale.linear().domain([0, d3.max(data)]).range([h, 0]);
  
  var line = d3.svg.line()
    .x(function(d,i) { 
        return x(d.t); 
    })
    .y(function(d) { 
        return y(d.value); 
    });

  graph.append("svg:path").attr("class","function").attr("d", line(data));
  graph.selectAll("circle").data(chart.data).enter()
    .append("circle")
    .attr("cx", function(d){ return x(d.t); })
    .attr("cy", function(d){ return y(d.value); })
    .attr("r", "5");
}

function createCharts(charts){
  _.map(charts, function(chart, key){
    var graph = basePlot(chart,key);
    if(chart.type == "signal"){
      plot(chart,key, graph);
    } else if(chart.type == "spectrum"){
      plotSpectrum(chart,key,graph);
    }
  });
}

createCharts(charts);
