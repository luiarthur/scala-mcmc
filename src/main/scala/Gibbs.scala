package mcmc

import java.util.Calendar
import scala.annotation.tailrec

trait Gibbs extends MCMC {
  // State
  type State

  // Substate to monitor. By default, it is the entire state.
  // thin period is 1 by default
  type Substate1
  def deepcopy1(s:State): Substate1
  val thin1:Int = 1


  // Optional substates to monitor.
  // thin period is 0 by default. i.e. no thinning.
  type Substate2
  def deepcopy2(s:State): Option[Substate2] = None
  val thin2:Int = 0

  type Substate3
  def deepcopy3(s:State): Option[Substate3] = None
  val thin3:Int = 0


  // Update function
  def update(s: State, iter:Int, out:Output): Unit

  type Output = (List[Substate1], List[Substate2], List[Substate3])

  def gibbs(s: State, niter:Int, nburn:Int,
            printProgress:Boolean=true): Output = {

    @tailrec
    def engine(_niter:Int, _nburn:Int,
               _out:Output=(List(), List(), List())): Output = {
      // Update the current state
      val iter = niter + nburn - (_niter + _nburn)
      update(s, _niter, _out)

      // Print progress
      if (printProgress && (_niter + _nburn) % 100 == 0) {
        println(s"Progress: ${_niter + _nburn} iterations remaining -- ${Calendar.getInstance.getTime}")
      }

      // Repeat if still burning
      if (_nburn > 0) {
        engine(_niter, _nburn - 1, _out)
      } else if (_niter > 0) {
        lazy val newOut1:List[Substate1] = {
          if (thin1 > 0 && _niter % thin1 == 0) {
            deepcopy1(s) :: _out._1
          } else { _out._1 }
        }

        lazy val newOut2:List[Substate2] = if (thin2 > 0 && _niter % thin2 == 0) {
          deepcopy2(s) match {
            case Some(x) => x :: _out._2
            case None => _out._2
          }
        } else { _out._2 }

        lazy val newOut3:List[Substate3] = if (thin3 > 0 && _niter % thin3 == 0) {
          deepcopy3(s) match {
            case Some(x) => x :: _out._3
            case None => _out._3
          }
        } else { _out._3 }

        val newOut = (newOut1, newOut2, newOut3)
        engine(_niter - 1, 0, newOut)
      } else {
        _out
      }
    }

    engine(niter, nburn)
  }
}

