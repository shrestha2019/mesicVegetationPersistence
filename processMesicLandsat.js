

// Applies scaling factors.
function applyScaleFactors(image) { 
  var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2);
  // var thermalBands = image.select('ST_B.*').multiply(0.00341802).add(149.0);
  return image.addBands(opticalBands, null, true)
              // .addBands(thermalBands, null, true);
}


exports.applyScaleFactors = applyScaleFactors; 

/* 
Cloud mask Sentinel 2 images  
*/

// MASKING CLOUDS
function maskL457sr(image) {
  // Bit 0 - Fill
  // Bit 1 - Dilated Cloud
  // Bit 2 - Unused
  // Bit 3 - Cloud
  // Bit 4 - Cloud Shadow
  var qaMask = image.select('QA_PIXEL').bitwiseAnd(parseInt('11111', 2)).eq(0);
  var saturationMask = image.select('QA_RADSAT').eq(0);

  // Apply the scaling factors to the appropriate bands.
  //var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2);
  //var thermalBand = image.select('ST_B6').multiply(0.00341802).add(149.0);

  // Replace the original bands with the scaled ones and apply the masks.
  return image
      .updateMask(qaMask)
      .updateMask(saturationMask);
}

exports.maskL457sr = maskL457sr; 



function etm2oliSR(image){
  /* Function to transform Landsat ETM+ SR values to OLI values
  using OLS regressions found in Roy et al 2016
  'Characterization of Landsat-7 to Landsat-8 reflective wavelength
  and normalized difference vegetation index continuity' RSE
  */
  
  //blue (~0.48um)
  var blue = image.expression('0.0003 + 0.8474 * B1', {
    'B1': image.select('SR_B1')});
  
   //green (~0.56um)
  var green = image.expression('0.0088 + 0.8483 * B2', {
    'B2': image.select('SR_B2')});
    
  //red (~0.66um)
  var red = image.expression('0.0061 + 0.9047 * B3', {
    'B3': image.select('SR_B3')});
  
  //nir (~0.85um)
  var nir = image.expression('0.0412 + 0.8462 * B4', {
    'B4': image.select('SR_B4')});
    
  //swir1 (~1.61um)
  var swir1 = image.expression('0.0254 + 0.8937 * B5', {
    'B5': image.select('SR_B5')});
  
  //swir2 (~2.21um)
  var swir2 = image.expression('0.0172 + 0.9071 * B7', {
    'B7': image.select('SR_B7')});
  
  //var bands = ee.List([blue, green, red, nir, swir1, swir2]);
  var bands = [blue, green, red, nir, swir1, swir2];
  
  // List of current band names
  var cur_names =  image.bandNames();
  
  // List of names of bands to be added
  var new_names = ee.List(['b', 'g', 'r','nir',  'swir1', 'swir2']);
  
  // flatten list to make one large list rather than list within the list
  var names = cur_names.add(new_names).flatten();
  
  return image.addBands(bands).rename(names);}





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
      'N': image.select('nir'),
      'R': image.select('r'),
      'G': image.select('g')
    }).rename('MCARI2');//Modified Chlorophyll Absorption in Reflectance Index 2
    
    var VSDI = image.expression(
      '1-(((SWIR)-Blue) + (Red -Blue))',
      // '1-((NIR-SWIR) + (NIR-Red))',
      {
      'SWIR': image.select('swir1'), 
      'Blue': image.select('b'), 
      'Red': image.select('r'), 
      'NIR': image.select('nir'), 
      'SWIR2': image.select('swir2'), 
      }).rename('VSDI')
  return image
    .addBands(ee.Image.constant(1))
    .addBands(ee.Image(years).rename('t')).float()
    .addBands(ee.Image(MCARI2))
    .addBands(ee.Image(VSDI))
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

function classify (image) {
  var mesics_MVsdi = image.expression(
          '((MCARI2 > 0.2 ) && (moist_VSDI > 0.4)) ? 1' +
          ': 0',
          
        {
          'moist_VSDI': image.select('moist_VSDI'),
          'MCARI2': image.select('MCARI2'),
          
        }).rename('mesics_MVsdi')


  return image
    .addBands(ee.Image(mesics_MVsdi))
}
 
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


var processL578 = function(startYear, endYear, startMonth, endMonth, geometry) {
  var LSbands=['b', 'g', 'r','nir',  'swir1', 'swir2']
  var imgsLS5 = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2")
        .filterBounds(geometry)
        .filter(ee.Filter.calendarRange(5,10,"month"))
        .filter(ee.Filter.calendarRange(1984, 2011, "year"))
        .filter(ee.Filter.lte('CLOUD_COVER', 20))
        .map(applyScaleFactors)
  var imgsLS7 = ee.ImageCollection("LANDSAT/LE07/C02/T1_L2")
        .filterBounds(geometry)
        .filter(ee.Filter.calendarRange(5,10,"month"))
        .filter(ee.Filter.calendarRange(2004, 2014, "year"))
        .filter(ee.Filter.lte('CLOUD_COVER', 20))
        .map(applyScaleFactors)
  var imgsLS = imgsLS5.merge(imgsLS7)
  var imgsLS57 = imgsLS.map(maskL457sr)
  // Apply transformation 
  imgsLS57 = imgsLS57.map(etm2oliSR).map(function(img){
    return img.select(['SR_B1','SR_B2','SR_B3','SR_B4','SR_B5','SR_B7', 'QA_PIXEL'])
    .rename('b', 'g', 'r','nir',  'swir1', 'swir2', 'QA_PIXEL')
    
  })
  // Landsat8 for 2013 - 2016
  var imgsLS8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
        .filterBounds(geometry)
        .filter(ee.Filter.calendarRange(5,10,"month"))
        .filter(ee.Filter.calendarRange(2013, 2022, "year"))
        .filter(ee.Filter.lte('CLOUD_COVER', 20))
        .map(applyScaleFactors)

//select and rename the bands to match images from 5 and 7
  imgsLS8 = imgsLS8.map(function(img){
  return img.select(['SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7', 'QA_PIXEL'])
          .rename('b', 'g', 'r','nir',  'swir1', 'swir2', 'QA_PIXEL')
})

// SR cloud AND cloud shadow mask
  var add_FmaskLS8 = function(image) {
    var msk = image.select('QA_PIXEL');
    msk = msk
    //Cloud
    .neq(21826).and(msk.neq(21890)).and(msk.neq(22080)).and(msk.neq(22144)).and(msk.neq(22280))
    //Shadow
    .and(msk.neq(23888)).and(msk.neq(23952)).and(msk.neq(24088)).and(msk.neq(24216)).and(msk.neq(24344)).and(msk.neq(24472))
    //Cirrus
    .and(msk.neq(54596)).and(msk.neq(54852)).and(msk.neq(55052))
    // Cirrus shadow
    .and(msk.neq(56856)).and(msk.neq(56984)).and(msk.neq(57240))
  return image.mask(msk);
    
  };

//Apply cloud mask
  var imgsLS8 = imgsLS8.map(add_FmaskLS8)

//Merge w 5 and 7
  var imgsLS578 = imgsLS57.merge(imgsLS8)
   var imgsLS578F = imgsLS578.filterDate(startYear+'-'+startMonth, endYear+'-'+endMonth)
  // var imgsLS578M = imgsLS578.map(maskLandform)//Mask image collection by landform/valley bottom
  
  var L578colIndices = imgsLS578F.map(calc_indices)
  var detrendedL578 = detrend(L578colIndices)
  var L578coldiff =  calc_moist_vsdi(detrendedL578)    
  var classifyL578 =  L578coldiff.map(classify)
  // var mesic_area = Classify.map(mesic_area)
  return classifyL578
}

exports.processL578 = processL578; 