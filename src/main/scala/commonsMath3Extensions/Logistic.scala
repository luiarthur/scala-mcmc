package org.apache.commons.math3.distribution

import org.apache.commons.math3.random.{RandomDataGenerator, RandomGeneratorFactory, RandomGenerator}
import org.apache.commons.math3.distribution._
import org.apache.commons.math3.util.FastMath.{log, pow, exp}
import mcmc.SpecialFunctions.{sigmoid, logit, sech}

/*
  val rng = RandomGeneratorFactory.createRandomGenerator(new java.util.Random(0))
  val rdg = new RandomDataGenerator(rng)
*/

case class Logistic(mean:Double, scale:Double, rng:RandomGenerator=null) extends AbstractRealDistribution(rng) {

  def cumulativeProbability(x:Double):Double = sigmoid((x - mean) / scale)
  private def pdfStandardized(x:Double):Double = {
    0.25 * math.pow(sech(x / 2), 2)
  }

  def density(x:Double):Double = (mean, scale) match {
    case (0, 1) => pdfStandardized(x)
    case _ => pdfStandardized((x - mean) / scale) / scale
  }

  lazy val variance:Double = pow(scale * math.Pi, 2.0) / 3.0

  def getNumericalMean:Double = mean
  def getNumericalVariance:Double = variance

  def getSupportLowerBound:Double = Double.NegativeInfinity
  def getSupportUpperBound:Double = Double.PositiveInfinity

  def isSupportConnected:Boolean = true
  def isSupportLowerBoundInclusive:Boolean = false
  def isSupportUpperBoundInclusive:Boolean = false

  override def sample():Double = {
    lazy val u = rng.nextDouble
    lazy val stdLogisticSample = log(u) - log(1 - u)

    (mean, scale) match {
      case (0, 1) => stdLogisticSample
      case _ => mean + scale * stdLogisticSample
    }
  }
}

