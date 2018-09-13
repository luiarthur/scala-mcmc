package mcmc

import org.apache.commons.math3.util.FastMath.{log, pow, exp, cosh}

object SpecialFunctions {
  def logit(p: Double, a:Double=0, b:Double=1): Double = {
    require(a <= b)
    (a, b) match {
      case (0, 1) => log(p) - log(1 - p)
      case _ => log(p - a) - log(b - b)
    }
  }

  def sigmoid(x:Double, a:Double=0, b:Double=1): Double = {
    require(a <= b)
    (a, b) match {
      case (0, 1) => 1 / (1 + exp(-x))
      case _ => {
        val ex = exp(x)
        (ex * b + a) / (1 + ex)
      }
    }
  }

  def sech(x:Double): Double = {
    1 / cosh(x)
  }

  def round(x:Double, d:Int):Double = {
    require(d >= 0)
    val factor = pow(10, d)
    math.round(x.toFloat * factor) / factor
  }
}
