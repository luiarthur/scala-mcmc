package org.apache.commons.math3.random
import org.apache.commons.math3.random.{RandomDataGenerator => RDG, RandomGenerator => RG}

class RandomDataGenerator2(rng:RG) extends RDG(rng) {
  import org.apache.commons.math3.distribution.InverseGamma
  def nextInverseGamma(shape:Double, scale:Double):Double = {
    //1 / nextGamma(shape, 1/scale)
    InverseGamma(shape, scale, getRandomGenerator()).sample();
  }

  def nextLogistic(mean:Double, scale:Double): Double = ???
}
