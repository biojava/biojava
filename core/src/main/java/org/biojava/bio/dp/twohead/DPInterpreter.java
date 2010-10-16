/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */


package org.biojava.bio.dp.twohead;

import java.io.Serializable;
import org.biojava.bio.BioError;
import org.biojava.bio.dp.BackPointer;
import org.biojava.bio.dp.DP;
import org.biojava.bio.dp.EmissionState;
import org.biojava.bio.dp.IllegalTransitionException;
import org.biojava.bio.dp.ScoreType;
import org.biojava.bio.dp.State;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;

/**
 * @author Matthew Pocock
 * @author Thomas Down
 */
public class DPInterpreter implements CellCalculatorFactory, Serializable {
  private final DP dp;

  public DPInterpreter(DP dp) {
    this.dp = dp;
  }

  public CellCalculator forwards(ScoreType scoreType)
  throws IllegalSymbolException, IllegalAlphabetException, IllegalTransitionException {
    return new Forward(dp, scoreType);
  }

  public CellCalculator backwards(ScoreType scoreType)
  throws IllegalSymbolException, IllegalAlphabetException, IllegalTransitionException {
    return new Backward(dp, scoreType);
  }

  public CellCalculator viterbi(ScoreType scoreType, BackPointer terminal)
  throws IllegalSymbolException, IllegalAlphabetException, IllegalTransitionException {
    return new Viterbi(dp, scoreType, terminal);
  }


  private static class Forward implements CellCalculator {
    private final int[][] transitions;
    private final double[][] transitionScores;
    private final State[] states;
    private final State magicalState;

    public Forward(DP dp, ScoreType scoreType)
    throws IllegalSymbolException, IllegalAlphabetException, IllegalTransitionException {
      states = dp.getStates();

      transitions = dp.getForwardTransitions();
      transitionScores = dp.getForwardTransitionScores(scoreType);
      magicalState = dp.getModel().magicalState();
    }

    public void initialize(Cell [][] cells)
    throws
      IllegalSymbolException,
      IllegalAlphabetException,
      IllegalTransitionException
    {
      _calcCell(cells, true);
    }

    public void calcCell(Cell [][] cells)
    throws
      IllegalSymbolException,
      IllegalAlphabetException,
      IllegalTransitionException
    {
      _calcCell(cells, false);
    }

    public void _calcCell(Cell [][] cells, boolean initializationHack)
    throws
      IllegalSymbolException,
      IllegalAlphabetException,
      IllegalTransitionException
    {
      Cell curCell = cells[0][0];
      double[] curCol = curCell.scores;
      double[] emissions = curCell.emissions;
      //System.out.println("curCol = " + curCol);

     STATELOOP:
      for (int l = 0; l < states.length; ++l) {
        State curState = states[l];
        //System.out.println("State = " + states[l].getName());
        try {
          if(initializationHack && (curState instanceof EmissionState)) {
            if(curState == magicalState) {
              curCol[l] = 0.0;
            } else {
              curCol[l] = Double.NaN;
            }
            //System.out.println("Initialized state to " + curCol[l]);
            continue STATELOOP;
          }

          //System.out.println("Calculating weight");
          double[] sourceScores;
          double weight;
          if (! (curState instanceof EmissionState)) {
            weight = 0.0;
            sourceScores = curCol;
          } else {
            weight = emissions[l];
            //System.out.println("Weight " + emissions[l]);
            if(weight == Double.NEGATIVE_INFINITY || Double.isNaN(weight)) {
              curCol[l] = Double.NaN;
              continue STATELOOP;
            }
            int [] advance = ((EmissionState)curState).getAdvance();
            sourceScores = cells[advance[0]][advance[1]].scores;
            //System.out.println("Values from " + advance[0] + ", " + advance[1] + " " + sourceScores);
          }
          //System.out.println("weight = " + weight);

          int [] tr = transitions[l];
          double[] trs = transitionScores[l];

          /*for(int ci = 0; ci < tr.length; ci++) {
            System.out.println(
              "Source = " + states[tr[ci]].getName() +
              "\t= " + sourceScores[tr[ci]]
            );
          }*/

          // Calculate probabilities for states with transitions
          // here.

          // Find base for addition
          double constant = Double.NaN;
          double score = 0.0;
          for (
            int ci = 0;
            ci < tr.length;
            ci++
          ) {
            int trc = tr[ci];
            double trSc = sourceScores[trc];
            if(!Double.isNaN(trSc) && trSc != Double.NEGATIVE_INFINITY) {
              if(Double.isNaN(constant)) {
                constant = trSc;
              }
              double sk = trs[ci];
              if(!Double.isNaN(sk) && sk != Double.NEGATIVE_INFINITY) {
                score += Math.exp(trSc + sk - constant);
              }
            }
          }
          if(Double.isNaN(constant)) {
            curCol[l] = Double.NaN;
            //System.out.println("found no source");
          } else {
            curCol[l] = weight + Math.log(score) + constant;
          }
        } catch (Exception e) {
          throw new BioError(

            "Problem with state " + l + " -> " + states[l].getName(), e
          );
        } catch (BioError e) {
          throw new BioError(

            "Error  with state " + l + " -> " + states[l].getName(), e
          );
        }
      }
      /*for (int l = 0; l < states.length; ++l) {
        State curState = states[l];
        System.out.println(
          "State = " + states[l].getName() +
          "\t = " + curCol[l]
        );
      }*/
    }
  }

  private static class Backward implements CellCalculator {
    private final int[][] transitions;
    private final double[][] transitionScores;
    private final State[] states;
    private final State magicalState;

    public Backward(
      DP dp, ScoreType scoreType
    ) throws
      IllegalSymbolException,
      IllegalAlphabetException,
      IllegalTransitionException
    {
      states = dp.getStates();
      transitions = dp.getBackwardTransitions();
      transitionScores = dp.getBackwardTransitionScores(scoreType);
      magicalState = dp.getModel().magicalState();
    }

    public void initialize(Cell [][] cells)
    throws
      IllegalSymbolException,
      IllegalAlphabetException,
      IllegalTransitionException
    {
      _calcCell(cells, true);
    }

    public void calcCell(Cell [][] cells)
    throws
      IllegalSymbolException,
      IllegalAlphabetException,
      IllegalTransitionException
    {
      _calcCell(cells, false);
    }

    public void _calcCell(Cell [][] cells, boolean initializationHack)
    throws IllegalSymbolException, IllegalAlphabetException, IllegalTransitionException
    {
      Cell curCell = cells[0][0];
      double[] curCol = curCell.scores;

     STATELOOP:
      for (int l = states.length - 1; l >= 0; --l) {
        //System.out.println("State = " + states[l].getName());
        State curState = states[l];
        if(initializationHack && (curState instanceof EmissionState)) {
          if(curState == magicalState) {
            curCol[l] = 0.0;
          } else {
            curCol[l] = Double.NaN;
          }
          continue STATELOOP;
        }
        int [] tr = transitions[l];
        double[] trs = transitionScores[l];

        // Calculate probabilities for states with transitions
        // here.

            double[] sourceScores = new double[tr.length];
        for (int ci = 0; ci < tr.length; ++ci) {
          double weight;

          int destI = tr[ci];
          State destS = states[destI];
          Cell targetCell;
          if (destS instanceof EmissionState) {
            int [] advance = ((EmissionState)destS).getAdvance();
            targetCell = cells[advance[0]][advance[1]];
            weight = targetCell.emissions[destI];
            if(Double.isNaN(weight)) {
              curCol[l] = Double.NaN;
              continue STATELOOP;
            }
          } else {
            targetCell = curCell;
            weight = 0.0;
          }
          sourceScores[ci] = targetCell.scores[destI] + weight;
        }
        /*for(int ci = 0; ci < tr.length; ci++) {
          System.out.println(
            "Source = " + states[tr[ci]].getName() +
            "\t= " + sourceScores[ci]
          );
        }*/

        // Find base for addition
        double constant = Double.NaN;
        double score = 0.0;
        for(
          int ci = 0;
          ci < tr.length;
          ci++
        ) {
          double skc = sourceScores[ci];
          if(skc != Double.NEGATIVE_INFINITY && !Double.isNaN(skc)) {
            if(Double.isNaN(constant)) {
              constant = skc;
            }
            double sk = trs[ci];
            if(!Double.isNaN(sk) && sk != Double.NEGATIVE_INFINITY) {
              score += Math.exp(skc + sk - constant);
            }
          }
        }
        if(Double.isNaN(constant)) {
          curCol[l] = Double.NaN;
        } else {
          curCol[l] = Math.log(score) + constant;
          //System.out.println(curCol[l]);
        }
      }
      /*for (int l = 0; l < states.length; ++l) {
        State curState = states[l];
        System.out.println(
          "State = " + states[l].getName() +
          "\t = " + curCol[l]
        );
      }*/
    }
  }


  private class Viterbi implements CellCalculator {
    private final int[][] transitions;
    private final double[][] transitionScores;
    private final State[] states;
    private final BackPointer TERMINAL_BP;
    private final State magicalState;

    public Viterbi(DP dp, ScoreType scoreType, BackPointer terminal)
    throws
      IllegalSymbolException,
      IllegalAlphabetException,
      IllegalTransitionException
    {
      TERMINAL_BP = terminal;
      states = dp.getStates();
      transitions = dp.getForwardTransitions();
      transitionScores = dp.getForwardTransitionScores(scoreType);
      magicalState = dp.getModel().magicalState();
    }

    public void initialize(Cell[][] cells)
    throws
      IllegalSymbolException,
      IllegalAlphabetException,
      IllegalTransitionException
    {
      _calcCell(cells, true);
    }

    public void calcCell(Cell[][] cells)
    throws
      IllegalSymbolException,
      IllegalAlphabetException,
      IllegalTransitionException
    {
      _calcCell(cells, false);
    }

    public void _calcCell(Cell [][] cells, boolean initializationHack)
    throws
      IllegalSymbolException,
      IllegalAlphabetException,
      IllegalTransitionException
    {
      Cell curCell = cells[0][0];
      double[] curCol = curCell.scores;
      BackPointer[] curBPs = curCell.backPointers;
      double[] emissions = curCell.emissions;
      //System.out.println("Scores " + curCol);
     STATELOOP:
      for (int l = 0; l < states.length; ++l) {
        State curState = states[l];
            //System.out.println("State = " + l + "=" + states[l].getName());
        try {
          //System.out.println("trying initialization");
          if(initializationHack && (curState instanceof EmissionState)) {
            if(curState == magicalState) {
              curCol[l] = 0.0;
              curBPs[l] = TERMINAL_BP;
            } else {
              curCol[l] = Double.NaN;
              curBPs[l] = null;
            }
            //System.out.println("Initialized state to " + curCol[l]);
            continue STATELOOP;
          }

          double weight;
          double[] sourceScores;
          BackPointer[] oldBPs;
          if(! (curState instanceof EmissionState)) {
            weight = 0.0;
            sourceScores = curCol;
            oldBPs = curBPs;
          } else {
            weight = emissions[l];
            if(weight == Double.NEGATIVE_INFINITY || Double.isNaN(weight)) {
              curCol[l] = Double.NaN;
              curBPs[l] = null;
              continue STATELOOP;
            }
            int [] advance = ((EmissionState)curState).getAdvance();
            Cell oldCell = cells[advance[0]][advance[1]];
            sourceScores = oldCell.scores;
            oldBPs = oldCell.backPointers;
            //System.out.println("Looking back " + advance[0] + ", " + advance[1]);
          }
          //System.out.println("weight = " + weight);

          double score = Double.NEGATIVE_INFINITY;
          int [] tr = transitions[l];
          double[] trs = transitionScores[l];

          int bestK = -1; // index into states[l]
          for (int kc = 0; kc < tr.length; ++kc) {
            int k = tr[kc]; // actual state index
            double sk = sourceScores[k];

            /*System.out.println("kc is " + kc);
            System.out.println("with from " + k + "=" + states[k].getName());
            System.out.println("prevScore = " + sk);*/
            if (sk != Double.NEGATIVE_INFINITY && !Double.isNaN(sk)) {
              double t = trs[kc];
              //System.out.println("Transition score = " + t);
              double newScore = t + sk;
              if (newScore > score) {
                score = newScore;
                bestK = k;
                //System.out.println("New best source at " + kc + " is " + score);
              }
            }
          }
          if (bestK != -1) {
            curCol[l] = weight + score;
            /*System.out.println("Weight = " + weight);
            System.out.println("Score = " + score);
            System.out.println(
              "Creating " + states[bestK].getName() +
              " -> " + states[l].getName() +
              " (" + curCol[l] + ")"
            );*/
            try {
              State s = states[l];
              curBPs[l] = new BackPointer(
                s,
                oldBPs[bestK],
                curCol[l]
              );
            } catch (Throwable t) {
              throw new BioError(

                "Couldn't generate backpointer for " + states[l].getName() +
                " back to " + states[bestK].getName(), t
              );
            }
          } else {
            //System.out.println("No where to come from");;
            curBPs[l] = null;
            curCol[l] = Double.NaN;
          }
        } catch (Exception e) {
          throw new BioError(

            "Problem with state " + l + " -> " + states[l].getName(),e
          );
        } catch (BioError e) {
          throw new BioError(

            "Error  with state " + l + " -> " + states[l].getName(),e
          );
        }
      }
      /*System.out.println("backpointers:");
      for(int l = 0; l < states.length; l++) {
        System.out.print(states[l].getName() + "\t" + curCol[l] + "\t");
        BackPointer b = curBPs[l];
        if(b != null) {
          for(BackPointer bb = b; bb.back != bb; bb = bb.back) {
            System.out.print(bb.state.getName() + " -> ");
          }
          System.out.println("!");
        } else {
          System.out.print("\n");
        }
      }*/
      initializationHack = false;
    }
  }

  public static class Maker implements CellCalculatorFactoryMaker {
    public CellCalculatorFactory make(DP dp) {
      return new DPInterpreter(dp);
    }
  }
}
