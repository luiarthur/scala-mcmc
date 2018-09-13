import org.apache.commons.math3.distribution._
import org.apache.commons.math3.random._
import mcmc.util.timer
import mcmc.SpecialFunctions.round

import org.scalactic.TolerantNumerics
import org.scalatest.FunSuite

class TestDistributions extends TestUtil {
  val printDebug = false
  val eps = 1E-8
  implicit val doubleEq = TolerantNumerics.tolerantDoubleEquality(eps)

  //val rng = RandomGeneratorFactory.createRandomGenerator(java.util.concurrent.ThreadLocalRandom.current())
  val rng = RandomGeneratorFactory.createRandomGenerator(new java.util.Random(0))

  test("Inverse-Gamma") {
    val ig = InverseGamma(3,5,rng)
    val x = 4.0
    assert(ig.density(x) === 0.0699474601)
    assert(ig.cumulativeProbability(x) ===0.8684676654)
    assert(ig.mean === 2.5)
    assert(ig.variance === 6.25, 1E-5)
    assert(ig.inverseCumulativeProbability(.6) === 2.188110164362591)
    val samps = ig.sample(1E6.toInt) // takes less than a second
    assert(round(samps.sum / samps.size, 2) === ig.mean)
  }

  test("Logistic") {
    val l = Logistic(3,2,rng)
    val x = 4.0
    assert(l.density(x) === 0.11750185610079725)
    assert(l.cumulativeProbability(x) ===0.6224593312018546)
    assert(l.mean === 3.0)
    assert(l.variance === 13.159472534785811)
    assert(l.inverseCumulativeProbability(.6) === 3.8109302162163283)
    val samps = l.sample(1E6.toInt) // takes less than a second
    assert(round(samps.sum / samps.size, 1) === l.mean)
  }
}
