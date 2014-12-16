
// Copyright 2013 Peter Cook @prcweb prcweb.co.uk
var chart = {

	data: null,
  xScale: null,
  yScale: null,
  svgLine: null,
  colorScale: null,

  chartHeight: 500,
  lineWidth: 600,
  lineHeight: 150,

  bodyHeight: 4000, // Scroll height
  windowHeight: 0,
  scrollScale: null,

  menu: [
    {'label': 'Year', 'sortBy': 'year'},
    {'label': 'Maximum', 'sortBy': 'max'},
    {'label': 'Minimum', 'sortBy': 'min'},
    {'label': 'Mean', 'sortBy': 'mean'}
  ],

  uiState: {
    selectedIndex: 0,
    selectedDatum: null, // Automatically updated
    sortBy: 'year',
    sorting: false
  },

  sortFunction: {},
  openingTimer: null,


  translate: function(x, y) {return 'translate('+x+','+y+')';},

  init: function() {
  	var data = {};

  	var	dataset = _.map(this.data.meanTemp, function(data, k) {
      var yearAve = _.reduce(data, function(m, v) {return m + v;}, 0) / 12;
      var yearMax = _.max(data);
      var yearMin = _.min(data);
			return {year: +k, data: data, mean: yearAve, max: yearMax, min: yearMin};
		});
  	this.data.meanTemp = {data: dataset.reverse() , extent: [-2, 16]};


  	d3.select('body').style('height', this.bodyHeight + 'px');
  	this.windowHeight = $(window).height();
  	this.scrollScale = d3.scale.linear().domain([0, this.bodyHeight - this.windowHeight * 0.99]).range([0, 102]).clamp(true);

    this.sortFunction.year = function(a, b) {return d3.descending(a.year, b.year);}
    this.sortFunction.mean = function(a, b) {return d3.descending(a.mean, b.mean);}
    this.sortFunction.max = function(a, b) {return d3.descending(a.max, b.max);}
    this.sortFunction.min = function(a, b) {return d3.ascending(a.min, b.min);}

  	this.initChart();
  	this.initEvents();
    this.initMenu();
  },

  initMenu: function() {
    var that = this;
    d3.select('#menu')
      .text('Sort by: ')
      .selectAll('span')
      .data(this.menu)
      .enter()
      .append('span')
      .html(function(d, i) {
        var html = '<span class="button">' + d.label + '</span>';
        if(i < that.menu.length - 1) html += ' / ';
        return html;
      });

    d3.select('#menu')
      .selectAll('span.button')
      .classed('selected', function(d, i) {return i===0;})
      .on('click', function() {
        var d = d3.select(this.parentNode).datum();
        console.log(d, d.sortBy);

        d3.selectAll('#menu span.button')
          .classed('selected', false);

        d3.select(this)
          .classed('selected', true);

        that.updateSort(d.sortBy);
      });
  },

  updateVisibleYears: function() {
  	var that = chart; // Better way to do this?

    var index = that.uiState.selectedIndex;
    var years = d3.selectAll('#chart .years g.year');
    years.classed('hover', false);

    years
      .filter(function(d, i) {return i === index;})
			.classed('hover', true);

		d3.selectAll('.axes')
  		.attr('transform', that.translate(30 + index * that.perspectiveOffsetX, that.chartHeight + that.yScale(0) + -index * that.perspectiveOffsetY));

    years
			.style('opacity', function(d, i) {
        if(i < index) return 0;
        return that.colorScale(i);
			});

    var datum = years.filter(function(d, i) {return i === index;}).datum();
    that.uiState.selectedDatum = datum;

    that.updateInfo();
  },

  updateInfo: function() {
    var that = chart;
    var d = that.uiState.selectedDatum;
    var html = '<h2>' + d.year + '</h2>';
    html += _.has(that.data.info, d.year) ? that.data.info[d.year].text : '';
    html += '<p>Hottest month: ' + d.max + '&deg;C</p>';
    html += '<p>Coolest month: ' + d.min + '&deg;C</p>';
    html += '<p>Year average: ' + d.mean.toFixed(1) + '&deg;C</p>';

    d3.select('#info')
      .html(html);
  },

  handleScroll: function() {
  	var that = chart; // Better way to do this?
    if(that.uiState.sorting) return;
		var scroll = $(window).scrollTop();
		that.uiState.selectedIndex = Math.round(that.scrollScale(scroll));
		that.updateVisibleYears();
  },

  initEvents: function() {
  	var that = this;
		$(window).scroll(this.handleScroll);
    $(window).on('touchmove', this.handleScroll);
  },

  initChart: function() {
  	var that = this;

  	this.xScale = d3.scale.linear()
  		.domain([0, 11])
  		.range([0, this.lineWidth]);

  	this.yScale = d3.scale.linear()
  		.domain(this.data.meanTemp.extent)
  		.range([this.lineHeight, 0]);

  	this.colorScale = d3.scale.linear()
  		.domain([0, 102])
  		.range([1, 0.5]);

    this.svgLine = d3.svg.line()
      .interpolate('cardinal')
      .x(function(d, i) {return that.xScale(i);})
      .y(function(d) {return that.yScale(d);});

    // YEAR LINES
  	var years = d3.select('#chart svg')
  		.append('g')
  		.classed('years', true)
  		.attr('transform', this.translate(30, this.chartHeight))
      .selectAll('g.year')
      .data(this.data.meanTemp.data)
      .enter()
      .append('g')
      .attr('class', function(d, i) {return 'year-' + d.year;})
      .classed('year', true)
      .sort(this.sortFunction[this.uiState.sortBy])
      .attr('transform', function(d, i) {
        return that.translate(i * that.perspectiveOffsetX, -i * that.perspectiveOffsetY);
      })
      .style('opacity', function(d, i) {
        return that.colorScale(i);
      });

    // Add paths
    years
      .append('path')
      .attr('d', function(d, i) {
        return that.svgLine(d.data);
      });

    // Base and end lines
    years
      .append('line')
      .classed('base', true)
      .attr('x1', 0)
      .attr('y1', this.yScale(0))
      .attr('x2', this.xScale(11))
      .attr('y2', this.yScale(0));

    years
      .append('line')
      .classed('start', true)
      .attr('x1', 0)
      .attr('y1', this.yScale(0))
      .attr('x2', 0)
      .attr('y2', function(d) {return that.yScale(d.data[0]);});

    years
      .append('line')
      .classed('end', true)
      .attr('x1', this.xScale(11))
      .attr('y1', this.yScale(0))
      .attr('x2', this.xScale(11))
      .attr('y2', function(d) {return that.yScale(d.data[11]);});

    years
      .append('text')
      .classed('label', true)
      .attr('x', this.xScale(11) + 5)
      .attr('y', this.yScale(0))
      .text(function(d) {return d.year;})
      .each(function(d) {
        if(!_.has(that.data.info, d.year)) return;
        var color = that.data.info[d.year].class === 'hot' ? 'indianred' : 'steelblue';
        d3.select(this)
          .style('opacity', 1)
          .style('font-weight', '100')
          .style('fill', color);
      });

  	d3.select('#chart svg')
  		.append('g')
  		.classed('axes', true)
  		.attr('transform', this.translate(30, this.chartHeight));

  	this.renderAxes();
  },

  renderAxes: function() {
		var monthScale = d3.scale.ordinal()
	    .domain(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
	    .rangePoints([0, this.lineWidth]);

    var yAxis = d3.svg.axis()
      .scale(this.yScale)
      .orient('left')
      .tickValues([0, 2, 4, 6, 8, 10, 12, 14, 16]);

    d3.select('#chart .axes')
      .append('g')
      .classed('axis y', true)
      .attr('transform', this.translate(0, -this.yScale(0)))
      .call(yAxis);

    var xAxis = d3.svg.axis()
      .scale(monthScale)
      .orient('bottom');

    d3.select('#chart .axes')
      .append('g')
      .classed('axis x', true)
      .call(xAxis);
  },

  updateSort: function(sortBy) {
    var that = this;
    this.uiState.sortBy = sortBy;

    // Do the sort
    var years = d3.select('#chart .years')
      .selectAll('g.year')
      .sort(this.sortFunction[this.uiState.sortBy]);

    // Persist the chosen year: get the index of the chosen year
    d3.selectAll('#chart .years g.year')
      .each(function(d, i) {
        if(d.year === that.uiState.selectedDatum.year) index = i;
      });
    that.uiState.selectedIndex = index;

    // Transform the axes
    d3.selectAll('.axes')
      .transition()
      .duration(2000)
      .attr('transform', that.translate(25 + index * that.perspectiveOffsetX, that.chartHeight + that.yScale(0) + -index * that.perspectiveOffsetY))
      .each('end', function() {
        that.uiState.sorting = false;
      });

    // Transform the year paths
    d3.selectAll('#chart .years .year')
      .transition()
      .duration(2000)
      .attr('transform', function(d, i) {
        return that.translate(i * that.perspectiveOffsetX, -i * that.perspectiveOffsetY);
      })
      .style('opacity', function(d, i) {
        if(i < index) return 0;
        return that.colorScale(i);
      });

    // Reset scroll
    this.uiState.sorting = true;
    $(window).scrollTop(this.scrollScale.invert(index));
  }
}

d3.json('data/temp.json', function(temperatureData) {

	chart.data = temperatureData;
	chart.init();
	chart.updateVisibleYears();
});

