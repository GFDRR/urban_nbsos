import ee
import math

# ****************************************************************************************************************** #
# ******************* Module compiling the compositing and time series processing functionalities ****************** #
# ****************************************************************************************************************** #

# A functioned to generate a harmonized time series from a Sentinel-2 Image Collection.
# Harmonized means that the generated temporal aggregates are equally spaced in time,
# i.e. by the number of days specified by the "agg_interval" argument
def harmonizedTS(masked_collection, band_list, time_intervals, options):

  band_name = options.band_name or 'NDVI'
  agg_type = options.agg_type or 'median'

  # a wrapper function for stacking the generated Sentinel-2 temporal aggregates
  def _stackBands(time_interval, stack):
    outputs = exports.aggregateStack(masked_collection, band_list, time_interval,
                                         {'agg_type': agg_type, 'band_name': band_name})

    return ee.List(stack).add(ee.Image(outputs))

  # Initialize the list of Sentinel_2 images.
  stack = ee.List([])

  # Populate the stack of sentinel-2 images with the specified bands/VIs
  agg_stack = ee.List(time_intervals).iterate(_stackBands, stack)

  return ee.ImageCollection(ee.List(agg_stack)).sort('system:time_start')

#Define function to extract time intervals to use to generate the temporal composites from Sentinel collections
def extractTimeRanges(start, end, agg_interval):
    """
    Extract the time range data from the received time range and aggregation interval e.g.,
    'input time interval': time_interval = ['2019-01-01', '2020-01-01'], 'agg_interval': 60 days
    'generate the following time intervals':
    time_range = [("2019-01-01T00:'00':00Z", "2019-03-01T00:'00':00Z"),
                ("2019-03-01T00:'00':00Z", "2019-05-01T00:'00':00Z"),
                ("2019-05-01T00:'00':00Z", "2019-07-01T00:'00':00Z"),
                ("2019-07-01T00:'00':00Z", "2019-09-01T00:'00':00Z"),
                ("2019-09-01T00:'00':00Z", "2019-11-01T00:'00':00Z"),
                ("2019-11-01T00:'00':00Z", "2020-01-01T00:'00':00Z")
    """

    start_date = ee.Date(start)
    end_date = ee.Date(end)

    # Number of intervals in the given "time_range" based on the specified "agg_interval" period
    interval_no = ee.Date(end).difference(ee.Date(start), 'day').divide(agg_interval).round()
    month_check = ee.Algorithms.If(ee.Number(30.4375 / agg_interval).round().gt(0),
                                       ee.Number(30.4375 / agg_interval).round(),
                                       ee.Number(1)); # The number of aggregation intervals within a month

    # Compute the relative date delta (in months) to add to each preceding period to compute the new one
    rel_delta = ee.Number(end_date.difference(start_date, 'day')) \
                    .divide(ee.Number(30.4375).multiply(interval_no)).ceil();

    # Compute the first time interval end date by adding the relative date delta (in months) to the start date
    end_date = start_date.advance(start_date.advance(rel_delta, 'month') \
                                  .difference(start_date, 'day') \
                                  .divide(month_check), 'day') \
                                  .advance(-1, 'second')

    time_intervals = ee.List([ee.List([start_date, end_date])])

    def iter_function(x, previous):
        start_date = ee.Date(ee.List(ee.List(previous).reverse().get(0)).get(1)).advance(1, 'second')
        end_date = start_date.advance(start_date.advance(rel_delta, 'month') \
                     .difference(start_date, 'day') \
                     .divide(month_check), 'day') \
                     .advance(-1, 'second')
        return ee.List(previous).add(ee.List([start_date, end_date]))

    time_intervals = ee.List(ee.List.sequence(1, interval_no.subtract(1)).iterate(iter_function, time_intervals))

    return time_intervals

# Define function to generate the temporally-aggregated image for a given "time_interval"
def aggregateStack(masked_collection, band_list, time_interval, options):

    band_name = options.band_name or 'NDVI'
    agg_type = options.agg_type or 'median'

    time_interval = ee.List(time_interval)
    agg_interval =  ee.Date(time_interval.get(1)).difference(time_interval.get(0), 'day')

    # Set the centre of the time interval as the system:time_start date
    timestamp = {'system:time_start': ee.Date(time_interval.get(0)) \
                                          .advance(ee.Number(agg_interval.divide(2)).ceil(), 'day') \
                                          .millis()}

    def iter_function(band, stack):
        return ee.Image(stack).addBands(ee.Image(0).mask())

    if agg_type == 'geomedian':
        # Reduces the time interval using the geomedian, as it performs better for auto-correlated variables,
        #  i.e. spectral bands. See: Roberts, D., Mueller, N., & McIntyre, A. (2017).
        # High-dimensional pixel composites from earth observation time series.
        # IEEE Transactions on Geoscience and Remote Sensing, 55(11), 6254-6264.
        # A condition is provided in case the time interval does not contain any cloud-free images,
        # which likely happens in the sub-tropics when a short aggregation interval is used.

        agg_image = ee.Algorithms.If(
          masked_collection.filterDate(time_interval.get(0), time_interval.get(1)).size().gt(0),
                                       masked_collection.filterDate(time_interval.get(0), time_interval.get(1)) \
                                       .select(band_list) \
                                       .reduce(ee.Reducer.geometricMedian(band_list.length), 4) \
                                       .rename(band_list) \
                                       .set(timestamp),
                                       ee.Image(ee.List(band_list.slice(1)) \
                                       .iterate(iter_function(),
                                                ee.Image(0).mask())).rename(band_list).set(timestamp))
    elif agg_type == 'median' :
        agg_image = ee.Algorithms.If(
          masked_collection.filterDate(time_interval.get(0), time_interval.get(1)).size().gt(0),
                                       masked_collection.filterDate(time_interval.get(0), time_interval.get(1)) \
                                       .select(band_list) \
                                       .reduce(ee.Reducer.median(), 4) \
                                       .rename(band_list) \
                                       .set(timestamp),
                                       ee.Image(ee.List(band_list.slice(1)) \
                                       .iterate(iter_function(),
                                                ee.Image(0).mask())).rename(band_list).set(timestamp))
    elif agg_type == 'max_band':
        agg_image = ee.Algorithms.If(
          masked_collection.filterDate(time_interval.get(0), time_interval.get(1)).size().gt(0),
                                       masked_collection.filterDate(time_interval.get(0), time_interval.get(1)) \
                                       .select(band_list) \
                                       .qualityMosaic(band_name) \
                                       .set(timestamp),
                                       ee.Image(ee.List(band_list.slice(1)) \
                                       .iterate(iter_function(),
                                                ee.Image(0).mask())).rename(band_list).set(timestamp))
    elif agg_type == 'sum':
        agg_image = ee.Algorithms.If(
          masked_collection.filterDate(time_interval.get(0), time_interval.get(1)).size().gt(0),
                                       masked_collection.filterDate(time_interval.get(0), time_interval.get(1)) \
                                       .select(band_list) \
                                       .reduce(ee.Reducer.sum(), 4) \
                                       .rename(band_list) \
                                       .set(timestamp),
                                       ee.Image(ee.List(band_list.slice(1)) \
                                                .iterate(iter_function(),
                                                         ee.Image(0).mask())).rename(band_list).set(timestamp))
    elif agg_type == 'count':
        agg_image = ee.Algorithms.If(
          masked_collection.filterDate(time_interval.get(0), time_interval.get(1)).size().gt(0),
                                       masked_collection.filterDate(time_interval.get(0), time_interval.get(1)) \
                                       .select(band_list) \
                                       .reduce(ee.Reducer.count(), 4) \
                                       .rename(band_list) \
                                       .set(timestamp),
                                       ee.Image(ee.List(band_list.slice(1)) \
                                                .iterate(iter_function(),
                                                         ee.Image(0).mask())).rename(band_list).set(timestamp))
    elif agg_type == 'mean':
              agg_image = ee.Algorithms.If(
                masked_collection.filterDate(time_interval.get(0), time_interval.get(1)).size().gt(0),
                                             masked_collection.filterDate(time_interval.get(0), time_interval.get(1)) \
                                             .select(band_list) \
                                             .reduce(ee.Reducer.mean(), 4) \
                                             .rename(band_list) \
                                             .set(timestamp),
                                             ee.Image(ee.List(band_list.slice(1)) \
                                                      .iterate(iter_function(),
                                                               ee.Image(0).mask())).rename(band_list).set(timestamp))


    return agg_image
