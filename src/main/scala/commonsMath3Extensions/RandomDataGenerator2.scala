package org.apache.commons.math3.random

import org.apache.commons.math3.random.{RandomDataGenerator => RDG, RandomGenerator => RG}
import org.apache.commons.math3.distribution.{InverseGamma, Logistic}

class RandomDataGenerator2(rng:RG) extends RDG(rng) {
  def nextInverseGamma(shape:Double, scale:Double):Double = {
    InverseGamma(shape, scale, getRandomGenerator()).sample()
  }

  def nextLogistic(mean:Double, scale:Double): Double = {
    Logistic(mean, scale, getRandomGenerator()).sample()
  }
}
