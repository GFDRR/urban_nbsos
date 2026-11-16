import ee

# ****************************************************************************************************************** #
# *********** Computing the different parameters of the Revised Universal Soil Loss Equation (RUSLE) *************** #
# ****************************************************************************************************************** #

# Computes the Soil Erodibility Factor K from: Renard, K., Foster, G., Weesies, G., McCool, D. & Yoder, D.
# Predicting Soil Erosion by Water: a Guide to Conservation Planning with
# the Revised Universal Soil Loss Equation (RUSLE) (USDA-ARS, Washington, 1997).


"""
Code translated to python available from SOILWATCH
https://github.com/SoilWatch/soil-erosion-watch/blob/master/RUSLE_factors.js

JS CODE WRTITTEN BY WILLIAM OUELETTE
"""
def factorK(africa=True):
    """
    Import soil covariate datasets from SoilGrids available under EUPL 1.2 licence:
    de Sousa, L., Poggio, L., Batjes, N.H., Heuvelink, G.B.M., Kempen, B., Ribeiro, E., Rossiter, D.
    SoilGrids 2.0: producing quality-assessed soil information for the globe. Under submission to SOIL
    """
    def _lowResInputs():
        clay = ee.Image("projects/soilgrids-isric/clay_mean")
        clay = clay.select('clay_0-5cm_mean') \
             .add(clay.select('clay_5-15cm_mean')) \
             .add(clay.select('clay_15-30cm_mean')) \
             .divide(3) \
             .divide(10) \
             .rename('clay')

        sand = ee.Image("projects/soilgrids-isric/sand_mean")
        sand = sand.select('sand_0-5cm_mean') \
             .add(sand.select('sand_5-15cm_mean')) \
             .add(sand.select('sand_15-30cm_mean')) \
             .divide(3) \
             .divide(10) \
             .clamp(0, 20) \
             .rename('sand')

        silt = ee.Image("projects/soilgrids-isric/silt_mean")
        silt = silt.select('silt_0-5cm_mean') \
             .add(silt.select('silt_5-15cm_mean')) \
             .add(silt.select('silt_15-30cm_mean')) \
             .divide(3) \
             .divide(10) \
             .clamp(0, 70) \
             .rename('silt')

        OM = ee.Image("projects/soilgrids-isric/soc_mean")
        OM = OM.select('soc_0-5cm_mean') \
           .add(OM.select('soc_5-15cm_mean')) \
           .add(OM.select('soc_15-30cm_mean')) \
           .divide(3) \
           .divide(100) \
           .multiply(1.72) \
           .clamp(0, 4) \
           .rename('OM')

        bulk_density = ee.Image("projects/soilgrids-isric/bdod_mean")
        bulk_density = bulk_density.select('bdod_0-5cm_mean') \
                     .add(bulk_density.select('bdod_5-15cm_mean')) \
                     .add(bulk_density.select('bdod_15-30cm_mean')) \
                     .divide(3) \
                     .divide(100) \
                     .rename('bulk_density')

        return clay.addBands(sand).addBands(silt).addBands(OM).addBands(bulk_density)

    # Import soil covariate datasets from the Africa Soil and Agronomy Data Cube, under Creative Commons License:
    # Hengl, T., Miller, M. A., Krizan, J., Shepherd, K. D., Sila, A., Kilibarda, M., ... & Crouch, J. (2021).
    # African soil properties and nutrients mapped at 30 m spatial resolution using two-scale ensemble machine learning
    # Scientific reports, 11(1), 1-18.
    def _highResInputs():
        clay = ee.Image("projects/sat-io/open-datasets/iSDAsoil_Africa_30m/clay_content") \
                 .select('b1') \
                 .rename('clay')

        sand = ee.Image("projects/sat-io/open-datasets/iSDAsoil_Africa_30m/sand_content") \
                 .select('b1') \
                 .clamp(0, 20) \
                 .rename('sand')

        silt = ee.Image("projects/sat-io/open-datasets/iSDAsoil_Africa_30m/silt_content") \
                 .select('b1') \
                 .clamp(0, 70) \
                 .rename('silt')

        OM = ee.Image("projects/sat-io/open-datasets/iSDAsoil_Africa_30m/carbon_organic") \
               .select('b1') \
               .divide(10) \
               .multiply(1.72) \
               .clamp(0, 4) \
               .rename('OM')

        bulk_density = ee.Image("projects/sat-io/open-datasets/iSDAsoil_Africa_30m/bulk_density") \
                         .select('b1') \
                         .rename('bulk_density')

        return clay.addBands(sand).addBands(silt).addBands(OM).addBands(bulk_density).unmask(0)

    if africa:
        soil_covariates = ee.Image(_highResInputs())
    else:
        soil_covariates = ee.Image(_lowResInputs())

    clay = soil_covariates.select('clay')
    sand = soil_covariates.select('sand')
    silt = soil_covariates.select('silt')
    OM = soil_covariates.select('OM')
    bulk_density = soil_covariates.select('bulk_density')

    # Packing density equation and texture classes taken from: R.J. Jones, G. Spoor, A. Thomasson
    # Vulnerability of subsoils in Europe to compaction: a preliminary analysis,
    # Exp. Impact Prev. Subsoil Compact. Eur. Union, 73 (2003), pp. 131-143
    pack_density = bulk_density.add(clay.multiply(0.009))

    texture_class = clay.lt(18).And(sand.gt(65))\
                        .add((clay.lt(35).And(sand.gt(15))).Or(clay.gt(18).And(sand.gt(65))).multiply(2)) \
                        .add(clay.lt(35).And(sand.lt(15)).multiply(3)) \
                        .add(clay.gte(35).And(clay.lte(60)).multiply(4)) \
                        .add(clay.gt(60).multiply(5))

    # Adding an organic texture class based on an Organic Matter content is 20%, and clay content is 0%,
    # or when Organic Matter content is 30%, and clay content is 50%. The values in between are linearly interpolated.
    # Organic soil definition taken from: Huang, P. T., Patel, M., Santagata, M. C., & Bobet, A. (2009).
    # Classification of organic soils.
    organic_soil_test = clay.updateMask(OM.gte(20)).interpolate([0, 50], [20, 30], 'clamp')
    texture_class = texture_class \
                    .add(texture_class.eq(0).And(OM.gte(organic_soil_test)).unmask(0).multiply(9));

    # Texture and packing density used as proxy for the calculation of the soil structure class, as done in
    # Panagos, P., Meusburger, K., Ballabio, C., Borrelli, P., & Alewell, C. (2014).
    # Soil erodibility in Europe: A high-resolution dataset based on LUCAS.
    # Science of the total environment, 479, 189-200.
    structure_class = texture_class.eq(1).And(pack_density.lt(1.40)).multiply(4) \
                          .add(texture_class.eq(1).And(pack_density.gte(1.40).And(pack_density.lte(1.75))).multiply(3)) \
                          .add(texture_class.eq(1).And(pack_density.gt(1.75)).multiply(2)) \
                          .add(texture_class.eq(2).And(pack_density.lt(1.40)).multiply(3)) \
                          .add(texture_class.eq(2).And(pack_density.gte(1.40).And(pack_density.lte(1.75))).multiply(2)) \
                          .add(texture_class.eq(2).And(pack_density.gt(1.75)).multiply(2)) \
                          .add(texture_class.eq(3).And(pack_density.lt(1.40)).multiply(2)) \
                          .add(texture_class.eq(3).And(pack_density.gte(1.40).And(pack_density.lte(1.75))).multiply(2)) \
                          .add(texture_class.eq(3).And(pack_density.gt(1.75)).multiply(1)) \
                          .add(texture_class.eq(4).And(pack_density.lt(1.40)).multiply(2)) \
                          .add(texture_class.eq(4).And(pack_density.gte(1.40).And(pack_density.lte(1.75))).multiply(1)) \
                          .add(texture_class.eq(4).And(pack_density.gt(1.75)).multiply(1)) \
                          .add(texture_class.eq(5).And(pack_density.lt(1.40)).multiply(2)) \
                          .add(texture_class.eq(5).And(pack_density.gte(1.40).And(pack_density.lte(1.75))).multiply(1)) \
                          .add(texture_class.eq(5).And(pack_density.gt(1.75)).multiply(1)) \
                          .add(texture_class.eq(9).And(pack_density.lt(1.40)).multiply(4)) \
                          .add(texture_class.eq(9).And(pack_density.gte(1.40).And(pack_density.lte(1.75))).multiply(3))

    # Global Hydrologic Soil Groups and Saturated Hydraulic Conductivity Ksat from:
    # Simons, G., Koster, R., & Droogers, P. (2020).
    # HiHydroSoil v2. 0-High Resolution Soil Maps of Global Hydraulic Properties.
    # Under MIT License, courtesy of: Copyright (c) 2021 Samapriya Roy
    hydrologic_group = ee.Image('projects/sat-io/open-datasets/HiHydroSoilv2_0/Hydrologic_Soil_Group_250m')
    Ksat = ee.Image('projects/sat-io/open-datasets/HiHydroSoilv2_0/ksat/Ksat_0-5cm_M_250m') \
               .add(ee.Image('projects/sat-io/open-datasets/HiHydroSoilv2_0/ksat/Ksat_5-15cm_M_250m')) \
               .add(ee.Image('projects/sat-io/open-datasets/HiHydroSoilv2_0/ksat/Ksat_15-30cm_M_250m')) \
               .divide(3) \
               .divide(1e4)

    cfvo = ee.Image('projects/soilgrids-isric/cfvo_mean')
    cfvo = cfvo.select('cfvo_0-5cm_mean') \
           .add(cfvo.select('cfvo_5-15cm_mean')) \
           .add(cfvo.select('cfvo_15-30cm_mean')) \
           .divide(3) \
           .divide(10000);

    # Adjusted Satured Hydraulic Conductivity from:
    # Brakensiek, D. L., Rawls, W. J., & Stephenson, G. R. (1986).
    # Determining the saturated hydraulic conductivity of a soil containing rock fragments.
    # Soil Science Society of America Journal, 50(3), 834-835.
    Ksat = Ksat.multiply(ee.Image(1).subtract(cfvo))

    # Mapping Global Hydrologic Soil Groups to their permeability class equivalents,
    # from high run-off potential (=~ low permeability), to low run-off potential (=~ high permeability).
    # Taken from: USDA (1983) National Soil Survey Handbook. No. 430, US Department of Agriculture, USDA, Washington DC
    permeability_class = Ksat.gt(146.304).And(hydrologic_group.eq(1))  \
                             .add(Ksat.lte(146.304).And(Ksat.gt(48.768)) \
                             .And(hydrologic_group.eq(1).Or(hydrologic_group.eq(14))).multiply(2)) \
                             .add(Ksat.lte(48.768).And(Ksat.gt(12.192)) \
                             .And(hydrologic_group.eq(2).Or(hydrologic_group.eq(24))).multiply(3)) \
                             .add(Ksat.lte(12.192).And(Ksat.gt(4.8768)) \
                             .And(hydrologic_group.eq(3)).multiply(4)) \
                             .add(Ksat.lte(4.8768).And(Ksat.gt(2.4384)) \
                             .And(hydrologic_group.eq(34)).multiply(5)) \
                             .add(Ksat.lte(2.4384) \
                             .And(hydrologic_group.eq(4)).multiply(6))

    M = sand.add(silt).multiply(ee.Image(100).subtract(clay)).rename('M') # Textural factor M

    # Soil Erodibility factor K final equation
    K = ee.Image(2.1).multiply(1e-4).multiply(M.pow(1.14)).multiply(ee.Image(12).subtract(OM)) \
            .add(ee.Image(3.25).multiply(structure_class.subtract(2))) \
            .add(ee.Image(2.5).multiply(permeability_class.subtract(3))) \
            .divide(100) \
            .multiply(0.1317)
    return K

# Retrieves the Global Rainfall Erosivity factor R taken from: Panagos, Panos, et al.
# "Global rainfall erosivity assessment based on high-temporal resolution rainfall records."
# Scientific reports 7.1 (2017): 1-12.
# The data stored in GEE assets is not exposed publicly due to strict terms of usage of ESDAC,
# whereby the data can under NO CIRCUMSTANCES be passed to third parties. It can be requested at the following link:
# https:#esdac.jrc.ec.europa.eu/content/global-rainfall-erosivity#tabs-0-description=1
def factorR():
    R = ee.Image('users/soilwatch/GlobalRainfallErosivity2010')
    return R

# Retrieves the Slope Length (L) and Slope Steepness (S) factors and combines them into a single LS factor.
def factorLS(slope_deg, slope_rad, slope_aspect):
    # Import upstream drainage area (referred to as contributing or accumulated area in soil erosion studies) from:
    # Yamazaki D., D. Ikeshima, J. Sosa, P.D. Bates, G.H. Allen, T.M. Pavelsky.
    # MERITHydro:A high-resolution global hydrography map based on latest topography datasets Water Resources Research,
    # vol.55, pp.5053-5073, 2019
    contrib_area = ee.Image("MERIT/Hydro/v1_0_1").select('upa') \
                       .multiply(1000000) \
                       .reproject('EPSG:4326', None, 30) \
                       .divide(9) \
                       .clamp(0,4000)
                                       # For instance: Zhang, Hongming, et al. An improved method for calculating
                                       # slope length (lambda) and the LS parameters of the Revised Universal Soil Loss
                                       # Equation for large watersheds. Geoderma 308 (2017):36-45

    # Calculate the slope direction, either 1 for 4-neighbours, or sqrt(2) for diagonal neighbours.
    slope_direction = slope_aspect.gt(337.5).And(slope_aspect.lte(360)) \
                                                .Or(slope_aspect.gte(0).And(slope_aspect.lte(22.5))).multiply(1) \
                          .add(slope_aspect.gt(22.5).And(slope_aspect.lte(77.5)).multiply(ee.Image(2).sqrt())) \
                          .add(slope_aspect.gt(77.5).And(slope_aspect.lte(112.5)).multiply(1)) \
                          .add(slope_aspect.gt(112.5).And(slope_aspect.lte(157.5)).multiply(ee.Image(2).sqrt())) \
                          .add(slope_aspect.gt(157.5).And(slope_aspect.lte(202.5)).multiply(1)) \
                          .add(slope_aspect.gt(202.5).And(slope_aspect.lte(247.5)).multiply(ee.Image(2).sqrt())) \
                          .add(slope_aspect.gt(247.5).And(slope_aspect.lte(292.5)).multiply(1)) \
                          .add(slope_aspect.gt(292.5).And(slope_aspect.lte(337.5)).multiply(ee.Image(2).sqrt()))

    # Rill to interrill erosion factor beta of the Slope length factor L
    beta = slope_rad.sin() \
               .divide(0.0896) \
               .divide(ee.Image(3).multiply(slope_rad.sin().pow(0.8).add(0.56)).abs())
    # Slope exponential factor m related to beta
    m = beta.divide(beta.add(1))

    # Calculation of the slope length (L) and slope steepness (S) factor following: Desmet, P.; Govers, G.
    # A GIS procedure for automatically calculating the ULSE LS factor on topographically complex landscape units.
    #J. Soil Water Conserv. 1996, 51, 427-433.
    L = ee.Image(2).multiply(contrib_area) \
            .divide(ee.Image(2).multiply(30).multiply(slope_direction).multiply(22.13)) \
            .pow(m) \
            .multiply(m.add(1))

    # Computing the Slope Steepness factor S from: Renard, K., Foster, G., Weesies, G., McCool, D. & Yoder, D.
    # Predicting Soil Erosion by Water: a Guide to Conservation Planning with
    # the Revised Universal Soil Loss Equation (RUSLE) (USDA-ARS, Washington, 1997).
    S = ee.Algorithms.If(slope_rad.tan().lt(0.09), # Convert slope in radian to steepness (in %)
                             slope_rad.sin().multiply(10.8).add(0.03), # Mild slope condition < 9 %
                             slope_rad.sin().multiply(16.8).subtract(0.5)); # Steep slope condition >= 9%

    return L.multiply(S)

# Retrieves the so-called Sustainability Factor S = 1 / (V*L) ,
# the inverse of the C and P factors of the original RUSLE, to draw contrast with the other erosion-inducing factors
def factorS(median_image, bs_freq, fcover_ts, sr_band_scale):

  sr_band_scale = sr_band_scale or 1

  # Landscape factor L, computed as per Karydas & Panagos, 2018, using Sobel 3x3 edge-detection kernel for convolution
  # This increases the importance of intra-field features such as trees, fences, hedgerows, etc,
  # thus providing protection against soil erosion.
  L = ee.Image(1).add(median_image.select('B8').convolve({'kernel': ee.Kernel.sobel()}).abs() \
          .add(median_image.select('B8').convolve({'kernel': ee.Kernel.sobel().rotate(1)})).abs() \
          .divide(sr_band_scale).sqrt())

  # The monthly vegetation factor V, defined as per Karydas & Panagos, 2018,
  # was computed for each month in the specified time interval and the integral was calculated over the time period.
  # Bare soil frequency data scaled and clamped to the range [5,8]
  # was used as proxy for the land use parameter, to simulate the range of conditions from
  # a degraded cropland (high bare soil frequency, 5) to a sustainably managed grassland (low bare soil frequency, 8).
  fcover_arr = fcover_ts.toArray()

  pw_mean = fcover_arr.arraySlice(0,1).add(fcover_arr.arraySlice(0,0,-1)).divide(2)
  fcover_integ = pw_mean.arrayReduce('mean',[0]).abs().toArray() \
                            .arraySlice(0, 0, 1).arrayProject([0]).arrayFlatten([['array']])

  V = ee.Image(1).subtract(bs_freq).multiply(10).clamp(5, 8).multiply(fcover_integ.divide(sr_band_scale)).exp()
  # The final Sustainability Factor S = 1 / (V*L)
  S = ee.Image(1).divide(L.multiply(V)).rename('sustainability_factor')

  return S


