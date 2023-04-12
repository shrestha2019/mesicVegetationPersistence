/* 
Cloud mask Sentinel 2 images  
*/
 
function maskS2clouds(image) {
  var qa = image.select('QA60');

  // Bits 10 and 11 are clouds and cirrus, respectively.
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;

  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
      .and(qa.bitwiseAnd(cirrusBitMask).eq(0));

  return image.updateMask(mask).divide(10000)
              .copyProperties(image, image.propertyNames())
  ;
}
exports.maskS2clouds = maskS2clouds; 

/* 
Calculate indices 
*/
function calc_indices(image){
    var date = ee.Date(image.get('system:time_start'));
    // var date = image.date();
    var years = date.difference(ee.Date('1970-01-01'), 'year');
    var NDVI = image.normalizedDifference(['B8', 'B4']).float().rename('NDVI')
    var MCARI2 = image.expression(
      '(1.5 * (2.5 * (N - R) - 1.3 * (N - G))) / ((((2.0 * N + 1) ** 2) - (6.0 * N - 5 * (R ** 0.5)) - 0.5) ** 0.5)',
      {
      'N': image.select('B8'),
      'R': image.select('B4'),
      'G': image.select('B3')
    }).rename('MCARI2');//Modified Chlorophyll Absorption in Reflectance Index 2
    
    var VSDI = image.expression(
      '1-(((SWIR2)-Blue) + (Red -Blue))',
      // '1-((NIR-SWIR) + (NIR-Red))',
      {
       'SWIR': image.select('B11'), 
       'Blue': image.select('B2'), 
       'Red': image.select('B4'), 
       'NIR': image.select('B8A'), 
       'SWIR2': image.select('B12'), 
       }).rename('VSDI')
  var NSDI2 = image.expression(
      '((SWIR-SWIR2)/(SWIR2))',
      // '1-((NIR-SWIR) + (NIR-Red))',
      {
       'SWIR': image.select('B11'), 
       'SWIR2': image.select('B12'), 
       }).rename('NSDI2')
  return image
    .addBands(ee.Image.constant(1))
    .addBands(ee.Image(years).rename('t')).float()
    .addBands(ee.Image(MCARI2))
    .addBands(ee.Image(VSDI))
    .addBands(ee.Image(NSDI2))
    .addBands(ee.Image(NDVI))
}


exports.calc_indices = calc_indices;


/* 
Detrend using timeseries images  
*/

function detrend (collection){
  return collection.map(function(image){
  var independents = ee.List(['constant', 't']);
  var dependent = ee.String('VSDI');
  var trend = collection.select(independents.add(dependent))
    .reduce(ee.Reducer.linearRegression(independents.length(), 1));
  var coefficients = trend.select('coefficients')
    .arrayProject([0])
    .arrayFlatten([independents]);

  return ee.Image(image.addBands(image.select(dependent).subtract(
          image.select(independents).multiply(coefficients).reduce('sum'))
          .rename(dependent)))
          .copyProperties(image, ['system:time_start']);
});
}

exports.detrend = detrend; 

/* 
Calculate changes in moisture using percentiles   
*/

function calc_moist_vsdi (collection){
  var per5 = collection.reduce(ee.Reducer.percentile([5]));
  var per95 = collection.reduce(ee.Reducer.percentile([95]));
  return collection.map(function(img){
    var moist_VSDI = img.expression(
      '(VSDI - VSDI_p5) /(VSDI_p95-VSDI_p5)', 
      {
          'VSDI': img.select('VSDI'),
          'VSDI_p5': per5.select('VSDI_p5'),
          'VSDI_p95': per95.select('VSDI_p95'),
        }).rename('moist_VSDI')
        return img
        .addBands(ee.Image(moist_VSDI))
  });
  
}

exports.calc_moist_vsdi = calc_moist_vsdi; 


/* 
Classify/threshold based on relative moisture and MCARI2
*/



function classify (collection, mcari, mci) {
  // var val1 = parseFloat(mcari)
  // var val2 = ee.Number(Number((mci)))
  return collection.map(function(image){
    var mesics_MVsdi = image.expression(
      '((MCARI2 > mcari ) && (moist_VSDI > mci)) ? 1' +
          ': 0',
          
        {
          'moist_VSDI': image.select('moist_VSDI'),
          'MCARI2': image.select('MCARI2'),
          'mcari':mcari,
          'mci':mci
          
        }).rename('mesics_MVsdi')
        return image
        .addBands(ee.Image(mesics_MVsdi))
    
  });
        
}



// function classify (image) {
//   var mesics_MVsdi = image.expression(
//           '((MCARI2 > 0.2 ) && (moist_VSDI > 0.4)) ? 1' +
//           ': 0',
          
//         {
//           'moist_VSDI': image.select('moist_VSDI'),
//           'MCARI2': image.select('MCARI2'),
          
//         }).rename('mesics_MVsdi')

//   var mesics_vsdiOnly = image.expression(
//           '(moist_VSDI >0.5) ? 1' +
//           // ":(NDMI > 0.4) ? 1" +
//           ': 0',
          
//         {
//           'moist_VSDI': image.select('moist_VSDI'),
//           'MCARI2': image.select('MCARI2'),
//           'NDMI': image.select('NDMI'),
          
//         }).rename('mesics_vsdiOnly')
        
//   var mesicsNDVI = image.expression(
//           '(NDVI >0.3) ? 1' +
//           // ":(NDMI > 0.4) ? 1" +
//           ': 0',
          
//         {
//           'NDVI': image.select('NDVI'),
          
//         }).rename('mesicsNDVI')
//   // var mesics_masked = mesics.multiply(lfmask).rename('mesics_masked'); 

//   return image
//     .addBands(ee.Image(mesics_MVsdi))
//     .addBands(ee.Image(mesics_vsdiOnly))
//     .addBands(ee.Image(mesicsNDVI))
// }
 
exports.classify = classify; 

// function mesic_area(image, aoi){
//     var image1 = ee.Image(1).mask(image.select('mesics_MVsdi'));
//     var date = ee.Date(image.get('system:time_start')).format('YYYY-MM-DD');
//     var area_pxa = image1.multiply(ee.Image.pixelArea()) 
//                     .reduceRegion(ee.Reducer.sum(),
//                     aoi,10,null,null,false,1e13)
//                     .get('constant');
//     var area = ee.Number(area_pxa).divide(1e6);                 
//     return image.set({'area': area, 'date':date
//     });
// }

function maskLandform(image){
  var landforms = ee.Image('CSP/ERGo/1_0/US/landforms').select('constant');
  // var landforms = dataset.select('constant');
  var lfmask = landforms.eq(24).or(landforms.eq(34)).or(landforms.eq(41)).or(landforms.eq(42))
              .or(landforms.eq(31)).or(landforms.eq(32)).or(landforms.eq(33))
  return image.updateMask(lfmask).copyProperties(image, image.propertyNames())
}


var processS2 = function(startYear, endYear, startMonth, endMonth, geometry, mcari, mci) {
  var s2coll = ee.ImageCollection('COPERNICUS/S2_HARMONIZED')
                  // .filterDate('2017-01-01', '2022-12-31')
                  .filterDate(startYear+'-'+startMonth, endYear+'-'+endMonth)
                  .filterBounds(geometry)
                  .filter(ee.Filter.calendarRange(5,10,'month'))
                  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',10))
                  // .map(maskS2clouds)
                  // .map(maskLandform) 
  var s2colNoClouds = s2coll.map(maskS2clouds)
  
  var s2colIndices = s2colNoClouds.map(calc_indices)
  var detrended = detrend(s2colIndices)
  var s2coldiff =  calc_moist_vsdi(detrended)    
  var Classify =  classify(s2coldiff, mcari, mci )
  // var Classify =  s2coldiff.map(classify)
  // var mesic_area = Classify.map(mesic_area)
  return Classify
}

exports.processS2 = processS2; 

// function meanPpt(img){
//   var mPpt = img.reduceRegion(ee.Reducer.mean(), 
//               geometry,4000,null,null,false,1e13)
//               .get('ppt')
//   var date = ee.Date(img.get('system:time_start'));
//   return img.set({'meanPpt': mPpt, 'date': date});
// }
// exports.meanPpt = meanPpt; 

// var processPpt = function(startYear, endYear, startMonth, endMonth, geometry) {
//   var datasetPpt = ee.ImageCollection("OREGONSTATE/PRISM/AN81d")
//                   .filterDate(startYear+'-'+startMonth, endYear+'-'+endMonth)
//                   .select('ppt')
//   var meanPpt = datasetPpt.map(meanPpt);
//     return meanPpt
  
// }

// exports.processPpt = processPpt; 

