import ee

# ****************************************************************************************************************** #
# ***************************** Computing the different steps of the BioDivMap Methodology ************************* #
# ****************************************************************************************************************** #

# Computes the Soil Erodibility Factor K from: Renard, K., Foster, G., Weesies, G., McCool, D. & Yoder, D.
# Predicting Soil Erosion by Water: a Guide to Conservation Planning with
# the Revised Universal Soil Loss Equation (RUSLE) (USDA-ARS, Washington, 1997).

def raoQ(region, options):

  subset_size = options.subset_size or 2000
  window_size = options.window_size or 4

  def wrap(img):

    subsets = img.select(0).neighborhoodToArray(ee.Kernel.square(window_size/2))
    N = ee.Number(window_size+1).pow(2)


    def iter_1(i, img):
        def iter_2(j, img):
            out2 = ee.Image(img).add(subsets.subtract(subsets.arraySlice(0, ee.Number(i).int(), ee.Number(i).add(1).int()) \
                                 .arraySlice(1, ee.Number(j).int(), ee.Number(j).add(1).int()) \
                                 .arrayProject([0]).arrayFlatten([['community']])).abs().multiply(ee.Image(1).divide(N)))
            return out2

        out1 = ee.Image(img).add(
            ee.List.sequence(
                ee.Number(i).add(1),
                window_size-1
            ).iterate(
                iter_2,
                subsets.subtract(
                    subsets.arraySlice(0, ee.Number(i).int(), ee.Number(i).add(1).int()) \
                        .arraySlice(1, 1, 2) \
                        .arrayProject([0]) \
                        .arrayFlatten([['community']])
                ).abs().multiply(ee.Image(1).divide(N))
            )
        )
        return out1

    rao_q = ee.Image(ee.List.sequence(1, window_size-2).iterate(iter_1, subsets.subtract(subsets.arraySlice(0, 0, 1) \
                               .arraySlice(1, 1, 2) \
                               .arrayProject([0]) \
                               .arrayFlatten([['community']])
                       ).abs() \
                       .multiply(ee.Image(1).divide(N))
    ))

    return img.select(0).multiply(100).toByte().entropy(ee.Kernel.square(window_size/2)).rename('shannon') \
           .addBands(rao_q.arrayReduce({'reducer':ee.Reducer.sum().unweighted(), 'axes':[0, 1]}) \
                          .arrayProject([0]).arrayFlatten([['rao']]).divide(2)) \
           .addBands(img.select(1)) \
           .copyProperties(img, ['system:time_start'])

  return wrap

