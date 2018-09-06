import org.scalatest.FunSuite

class TestAdaptiveMetropolis extends FunSuite with mcmc.MCMC {
  object rng extends distribution.RandomSeq(new scala.util.Random(10))

  import mcmc.TuningParam


  def arrayToString[T](arr:Array[T]): String = arr.mkString(", ")
  def mean(x:List[Double]):Double = x.sum / x.size
  def variance(x:List[Double]):Double = {
    val xbar = mean(x)
    val n = x.size
    x.map{ xi => (xi - xbar) * (xi-xbar) / (n-1) }.sum
  }
  def sd(x:List[Double]):Double = math.sqrt(variance(x))

  val printDebug = false

  rng.setSeed(0)

  test("Normal Model") {
    val muTrue = Array(1.0, 3.0)
    val sig2True = 0.5
    val nHalf = 300
    val y = muTrue.map{ m => 
      Array.tabulate(nHalf){ i => rng.nextGaussian(m, math.sqrt(sig2True)) }
    }
    val n = nHalf * 2
    val muPriorMean = 0.0
    val muPriorSd = 5.0

    case class Param(var mu:Array[Double], var sig2:Double)
    object Model extends mcmc.Gibbs {
      type State = Param
      type Substate1 = Param
      val stepSigMu = Array.fill(muTrue.size){ TuningParam(1.0) }
      def deepcopy1(s:State):Substate1 = Param(s.mu.clone, s.sig2 + 0.0)
      def updateMuj(s:State, j:Int, i:Int, out:Output): Unit = {
        def logFullCond(muj:Double):Double = {
          val ll = y(j).map{ yj => -math.pow(yj - muj, 2.0) / (2*s.sig2) }.sum
          val lp = -math.pow(muj - muPriorMean, 2.0) / (2*muPriorSd*muPriorSd)
          ll + lp
        }

        s.mu(j) = metropolisAdaptive(s.mu(j), logFullCond, stepSigMu(j), rng)
      }

      def updateMu(s:State, i:Int, out:Output): Unit = {
        (0 to 1).foreach{ j => updateMuj(s, j=j, i, out) }
      }

      def update(s:State, i:Int, out:Output) {
        // updateSig2(s, i, out) // TODO 
        updateMu(s, i, out)
      }
    }

    val state = Param(Array(0, 0), 1)
    val (niter, nburn) = (2000, 2000)

    val out = Model.gibbs(state, niter=niter, nburn=nburn, printProgress=false)
    val muPost = out._1.map{ s => s.mu }
    val muMean = muPost.transpose.map{ mean }
    val muSd = muPost.transpose.map{ sd }

    if (printDebug) {
      out._1.foreach{ s => println(arrayToString(s.mu)) }
      println(s"mu post mean: ${muMean}")
      println(s"mu truth: ${muTrue.toList}")
      println(muSd)
    }

    val eps = 0.05
    muMean.toList.zip(muTrue).foreach{ case (m, mtrue) =>
      assert(math.abs(m - mtrue) < eps)
    }

    Model.stepSigMu.indices.foreach{ j => 
      println(s"mu(${j}) Acceptance Rate: ${Model.stepSigMu(j).accRate}")
    }

    assert(true)
  }
}


