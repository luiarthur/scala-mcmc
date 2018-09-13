import org.scalatest.FunSuite

class TestAdaptiveMetropolis extends TestUtil with mcmc.MCMC {
  val printDebug = false
  import mcmc.TuningParam

  // Set random seed
  object rng extends distribution.RandomSeq(new scala.util.Random(10))
  rng.setSeed(0)
  //import distribution.{RandomPar => rng}

  test("Normal Model") {
    import distribution.continuous.{InverseGamma, Normal}

    val muTrue = Array(1.0, 3.0)
    val J = muTrue.size
    val sig2True = 0.5
    val nHalf = 1000
    val y = muTrue.map{ m => 
      Array.tabulate(nHalf){ i => rng.nextGaussian(m, math.sqrt(sig2True)) }
    }
    val n = nHalf * 2
    val muPriorMean = 0.0
    val muPriorSd = 5.0
    val sig2PriorA = 3.0
    val sig2PriorB = 2.0

    case class Param(var mu:Array[Double], var sig2:Double)
    object Model extends mcmc.Gibbs {
      type State = Param
      type Substate1 = Param
      val stepSigMu = Array.fill(muTrue.size){ TuningParam(1.0) }
      val stepSigSig2 = TuningParam(1.0)
      def deepcopy1(s:State):Substate1 = Param(s.mu.clone, s.sig2)
      def updateMuj(s:State, j:Int): Unit = {
        def logFullCond(muj:Double):Double = {
          val ll = y(j).map{ Normal(muj, s.sig2).lpdf }.sum
          val lp = Normal(muPriorMean, muPriorSd).lpdf(muj)
          ll + lp
        }

        s.mu(j) = metropolisAdaptive(s.mu(j), logFullCond, stepSigMu(j), rng)
      }

      def updateMu(s:State): Unit = {
        (0 to 1).foreach{ j => updateMuj(s, j=j) }
      }

      def updateSig2(s:State): Unit = {
        def lp(sig2:Double):Double = InverseGamma(sig2PriorA,sig2PriorB).lpdf(sig2)
        def ll(sig2:Double):Double = {
          (0 until J).map{ j =>
            y(j).map{ yij => Normal(s.mu(j), math.sqrt(sig2)).lpdf(yij) }.sum
          }.sum
        }
        s.sig2 = metLogAdaptive(s.sig2, ll, lp, stepSigSig2, rng)
      }

      def update(s:State, i:Int, out:Output) {
        // updateSig2(s, i, out) // TODO 
        updateMu(s)
        updateSig2(s)
      }
    }

    val state = Param(Array(0, 0), 1)
    val (niter, nburn) = (2000, 2000)

    val out = Model.gibbs(state, niter=niter, nburn=nburn, printProgress=false)
    val muPost = out._1.map{ s => s.mu }
    val muMean = muPost.transpose.map{ mean }
    val muSd = muPost.transpose.map{ sd }
    val sig2Post = out._1.map{ _.sig2 }
    val sig2Mean = mean(sig2Post)

    if (printDebug) {
      out._1.foreach{ s => println(arrayToString(s.mu)) }
      println(s"mu post mean: ${muMean}")
      println(s"mu truth: ${muTrue.toList}")

      println(s"sig2 post mean: ${sig2Mean}")
      println(s"sig2 truth: ${sig2True}")
      println(muSd)

      Model.stepSigMu.indices.foreach{ j => 
        println(s"mu(${j}) Acceptance Rate: ${Model.stepSigMu(j).accRate}")
      }

      println(s"sig2 Acceptance Rate: ${Model.stepSigSig2.accRate}")
    }

    val eps = 0.05
    muMean.toList.zip(muTrue).foreach{ case (m, mtrue) =>
      assertApprox(m, mtrue, eps)
    }
    assertApprox(sig2Mean, sig2True, eps)

    /* Rscala example
    val R = org.ddahl.rscala.RClient()
    R.sig2Post = sig2Post.toArray
    R eval """
    library(rcommon)
    pdf("src/test/output/plots.pdf")
    plotPost(sig2Post)
    dev.off()
    """ 
    */
  }
}


