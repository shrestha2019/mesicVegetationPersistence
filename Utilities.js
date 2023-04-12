// Find all available NAIP images for a geometry
function findNAIP(geometry) {
  var init_collection = ee.ImageCollection('USDA/NAIP/DOQQ') 
    .filterBounds(geometry)
    .filterDate('2002-01-01', '2022-12-31')
    // .filter(ee.Filter.listContains("system:band_names", "N"));

  var yearList = ee.List(init_collection.distinct(['system:time_start']).aggregate_array('system:time_start'));
  var init_years = yearList.map(function(y){
    return ee.Date(y).get('year');
  });

  // remove duplicates
  init_years = ee.Dictionary(init_years.reduce(ee.Reducer.frequencyHistogram())).keys();
  var years = init_years.map(function(x) {return ee.Number.parse(x)});

  // Available NAIP years with NIR band
  var NAIPAnnual= function(year){
    var start_date = ee.Date.fromYMD(year, 1, 1);
    var end_date = ee.Date.fromYMD(year, 12, 31);
    var collection = init_collection
      .filterDate(start_date, end_date);

    var time_start = ee.List(collection.aggregate_array('system:time_start')).sort().get(0);
    var time_end = ee.List(collection.aggregate_array('system:time_end')).sort().get(-1);
    var col_size = collection.size();
    var image = ee.Image(collection.mosaic()
    // .clip(geometry)
    );

    return image.set({'system:time_start': time_start, 'system:time_end': time_end, 'years': year});
  };

  var naip = ee.ImageCollection(years.map(NAIPAnnual));

  return naip;
}

exports.findNAIP = findNAIP; 
