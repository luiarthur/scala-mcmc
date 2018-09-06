package mcmc

import distribution.RandomGeneric

case class TuningParam[T](var value:T, var accCount:Int=0, var currIter:Int=1) {

  def update(accept:Boolean) {
    if (accept) accCount += 1
    currIter += 1
  }

  def accRate:Double = accCount.toDouble / currIter
}


trait MCMC {
  type RNG = RandomGeneric
  def metropolis(curr:Double, logFullCond: Double=>Double,
                 stepSig:Double, rdg:RNG): Double = {

    val cand = rdg.nextGaussian(curr, stepSig)
    val u = math.log(rdg.nextUniform(0,1))
    val p = logFullCond(cand) - logFullCond(curr)
    if (p > u) cand else curr
  }

  def metropolisVec(curr:Array[Double], logFullCond:Array[Double]=>Double, stepCovMat:Array[Array[Double]], rdg:RNG): Array[Double] = {
    val cand = rdg.nextMvNormal(curr, stepCovMat)
    val u = math.log(rdg.nextUniform(0,1))
    val p = logFullCond(cand) - logFullCond(curr)
    if (p > u) cand else curr
  }

  /* Adaptive metropolis (within Gibbs). See section 3 of the paper below:
   *   http://probability.ca/jeff/ftpdir/adaptex.pdf
   *
   * Another useful website:
   *   https://m-clark.github.io/docs/ld_mcmc/index_onepage.html
   */
  def metropolisAdaptive(curr:Double, logFullCond:Double=>Double,
                         stepSig:TuningParam[Double], rdg:RNG,
                         delta:Int=>Double=defaultDelta,
                         targetAcc:Double=.44):Double = {

    val iter = stepSig.currIter
    val factor = math.exp(delta(iter))

    if (stepSig.accRate > targetAcc) {
      stepSig.value *= factor
    } else {
      stepSig.value /= factor
    }

    val cand = rdg.nextGaussian(curr, stepSig.value)
    val u = math.log(rdg.nextUniform(0,1))
    val p = logFullCond(cand) - logFullCond(curr)
    val accept = p > u

    stepSig.update(accept)

    if (accept) cand else curr
  }

  def defaultDelta(n:Int):Double = {
    List(math.pow(n.toDouble,-0.5), 0.01).min
  }

  def sech(x:Double):Double = 1 / math.cosh(x)

  def sigmoid(x:Double, a:Double=0, b:Double=1): Double = {
    (a, b) match {
      case (0, 1) => 1 / (1 + math.exp(-x))
      case _ => {
        val ex = math.exp(x)
        (b * ex + a) / (1 + ex)
      }
    }
  }

  def logit(p:Double, a:Double=0, b:Double=1): Double = {
    math.log(p - a) - math.log(b - p)
  }

  // TODO: Needs testing.
  def logpdfLogX2Param(logX:Double, logpdfX:(Double,Double,Double)=>Double, params:Double, a:Double, b:Double): Double = {
    logpdfX(math.exp(logX), a, b) + logX
  }

  def logpdfLogX(logX:Double, logpdfX:Double=>Double): Double = {
    logpdfX(math.exp(logX)) + logX
  }
  /* R Test
  logpdfLogX = function(logX, logpdfX) logpdfX(exp(logX)) + logX
  logX = log(rgamma(1E6, 5, 3))
  xx = seq(-10, 3, l=1000)
  den = logpdfLogX(xx, function(x) dgamma(x, 5, 3, log=TRUE))
  plot(density(logX), col='black', lwd=3)
  lines(xx, exp(den), lty=3, col='grey', lwd=3)
   */

  def logpdfLogitX(logitX:Double, logpdfX:Double=>Double, a:Double=0, b:Double=1): Double = {
    lazy val x = sigmoid(logitX, a, b)
    lazy val logJacobian:Double = logpdfLogistic(logitX) + math.log(b - a)
    logpdfX(x) + logJacobian
  }
  /* R Test
  sigmoid = function(y, a=0, b=1) (b * exp(y) + a) / (1 + exp(y))
  logit = function(x, a=0, b=1) log(x - a) - log(b - x)
  logpdfLogitX = function(logitX, logpdfX, a=0, b=1) {
    logpdfX(sigmoid(logitX, a, b)) + dlogis(logitX, log=TRUE) + log(b-a)
  }
  logitX = logit(runif(1E6, 2, 6), 2, 6)
  xx = seq(-10, 10, l=1000)
  den = logpdfLogitX(xx, function(x) dunif(x, 2, 6, log=TRUE), 2, 6)
  plot(density(logitX), col='black', lwd=3)
  lines(xx, exp(den), lty=3, col='red', lwd=3)
   */

  
  def pdfLogistic(x:Double, loc:Double=0, scale:Double=1):Double = {
    (loc, scale) match {
      case (0, 1) => 0.25 * math.pow(sech(x / 2.0), 2.0)
      case _ => pdfLogistic((x - loc) / scale) / scale
    }
  }

  def logpdfLogistic(x:Double, loc:Double=0, scale:Double=1): Double = {
    (loc, scale) match {
      case (0, 1) => math.log(0.25) + 2.0 * math.log(sech(x / 2.0))
      case _ => logpdfLogistic((x - loc) / scale) - math.log(scale)
    }
  }

  // TODO: Test
  def metLogAdaptive(curr:Double, ll:Double=>Double, lp: Double=>Double,
                     stepSig:TuningParam[Double], rdg:RNG,
                     delta:Int=>Double=defaultDelta,
                     targetAcc:Double=.44):Double = {
    val currLogX = math.log(curr)

    def lfcLogX(logX:Double):Double = {
      val x = math.exp(logX)
      ll(x) + logpdfLogX(logX, lp)
    }

    val newLogX = metropolisAdaptive(currLogX, lfcLogX, stepSig, rdg,
                                     delta, targetAcc)
    val newX = math.exp(newLogX)
    newX
  }

  // TODO: Test
  def metLogitAdaptive(curr:Double, ll:Double=>Double, lp: Double=>Double,
                       a:Double=0, b:Double=1,
                       stepSig:TuningParam[Double], rdg:RNG,
                       delta:Int=>Double=defaultDelta,
                       targetAcc:Double=.44):Double = {
    val currLogitX = logit(curr, a, b)

    def lfcLogitX(logitX:Double):Double = {
      val x = math.exp(logitX)
      ll(x) + logpdfLogitX(logitX, lp, a, b)
    }

    val newLogitX = metropolisAdaptive(currLogitX, lfcLogitX, stepSig, rdg,
                                       delta, targetAcc)
    val newX = sigmoid(newLogitX, a, b)
    newX
  }
}
