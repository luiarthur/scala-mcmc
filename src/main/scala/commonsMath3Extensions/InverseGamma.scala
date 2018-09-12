package org.apache.commons.math3.distribution

import org.apache.commons.math3.random.{RandomDataGenerator, RandomGeneratorFactory, RandomGenerator}
import org.apache.commons.math3.distribution._

/*
  val rng = RandomGeneratorFactory.createRandomGenerator(new java.util.Random(0))
  val rdg = new RandomDataGenerator(rng)
*/

case class InverseGamma(val shape:Double, val scale:Double, rng:RandomGenerator=null) extends AbstractRealDistribution(rng) {
  import math.{log, pow, exp}
  import org.apache.commons.math3.special.Gamma.{logGamma, regularizedGammaQ}
  override def logDensity(x:Double) = if (x > 0) {
    shape * log(scale) - logGamma(shape) - (shape+1) * log(x) - scale / x
  } else {
    Double.NegativeInfinity
  }

  def density(x:Double) = exp(logDensity(x))

  def cumulativeProbability(x:Double) = if (x > 0) {
    regularizedGammaQ(shape, scale / x)
  } else 0

  lazy val mean:Double = if (shape > 1) scale / (shape - 1) else Double.PositiveInfinity
  lazy val variance:Double = if (shape > 2) pow(scale / (shape - 1), 2) / (shape - 2) else Double.PositiveInfinity
  lazy val getNumericalMean: Double = mean
  lazy val getNumericalVariance: Double = variance
  lazy val mode = scale / (shape + 1)
  lazy val getSupportLowerBound: Double = 0.0
  lazy val getSupportUpperBound: Double = Double.PositiveInfinity
  lazy val isSupportConnected:Boolean = true
  lazy val isSupportLowerBoundInclusive:Boolean = false
  lazy val isSupportUpperBoundInclusive:Boolean = false
  override def sample():Double = {
    // Note that Gamma use (shape, scale) parameterization
    1 / (new GammaDistribution(rng, shape, 1 / scale).sample)
  }
}

