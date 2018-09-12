package org.apache.commons.math3.distribution

import org.apache.commons.math3.random.{RandomDataGenerator, RandomGeneratorFactory, RandomGenerator}
import org.apache.commons.math3.distribution._

/*
  val rng = RandomGeneratorFactory.createRandomGenerator(new java.util.Random(0))
  val rdg = new RandomDataGenerator(rng)
*/

case class Logistic(mean:Double, scale:Double, rng:RandomGenerator=null) extends AbstractRealDistribution(rng) {
  import math.{log, pow, exp}
}

