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
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.alignment.Alignment;
import org.biojava.bio.alignment.SimpleAlignment;
import org.biojava.bio.dp.BackPointer;
import org.biojava.bio.dp.DP;
import org.biojava.bio.dp.DPMatrix;
import org.biojava.bio.dp.EmissionState;
import org.biojava.bio.dp.IllegalTransitionException;
import org.biojava.bio.dp.MarkovModel;
import org.biojava.bio.dp.ScoreType;
import org.biojava.bio.dp.SimpleStatePath;
import org.biojava.bio.dp.State;
import org.biojava.bio.dp.StatePath;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.DoubleAlphabet;
import org.biojava.bio.symbol.GappedSymbolList;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.SimpleGappedSymbolList;
import org.biojava.bio.symbol.SimpleSymbolList;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.SmallMap;

/**
 * Algorithms for dynamic programming (alignments) between pairs
 * of SymbolLists.
 * Based on a single-head DP implementation by Matt Pocock.
 *
 * @author Thomas Down
 * @author Matthew Pocock
 */

public class PairwiseDP extends DP implements Serializable {
  private final HashMap emissionCache;
  private final CellCalculatorFactory ccFactory;

  public PairwiseDP(MarkovModel mm, CellCalculatorFactoryMaker ccfm)
  throws
    IllegalSymbolException,
    IllegalTransitionException,
    BioException
  {
    super(mm);
    Alphabet alpha = mm.emissionAlphabet();
    emissionCache = new HashMap();
    emissionCache.put(ScoreType.PROBABILITY, new EmissionCache(
      alpha, getStates(), getDotStatesIndex(), ScoreType.PROBABILITY)
    );
    emissionCache.put(ScoreType.ODDS, new EmissionCache(
      alpha, getStates(), getDotStatesIndex(), ScoreType.ODDS)
    );
    emissionCache.put(ScoreType.NULL_MODEL, new EmissionCache(
      alpha, getStates(), getDotStatesIndex(), ScoreType.NULL_MODEL)
    );
    this.ccFactory = ccfm.make(this);
  }
  
  private EmissionCache getEmissionCache(ScoreType scoreType) {
    return (EmissionCache) emissionCache.get(scoreType);
  }

  //
  // BACKWARD
  //

  public void update() {
    super.update();
    // workaround for bug in vm
    if(emissionCache != null) {
      for(Iterator i = emissionCache.values().iterator(); i.hasNext(); ) {
        ((EmissionCache) i.next()).clear();
      }
    }
  }

  private Cell run(PairDPCursor cursor, CellCalculator cc)
      throws IllegalSymbolException, IllegalAlphabetException, IllegalTransitionException
  {
    Cell[][] cells = cursor.press();
    if(cursor.hasNext()) {
      cursor.next(cells);
      cc.initialize(cells);
    }
    while(cursor.hasNext()) {
      cursor.next(cells);
      cc.calcCell(cells);
    }
    return cells[0][0];
  }
  
  private double runFB(PairDPCursor cursor, CellCalculator cc) 
      throws IllegalSymbolException, IllegalAlphabetException, IllegalTransitionException
  {
    Cell cell = run(cursor, cc);
    
    // Terminate!
    State [] states = getStates();
    int l = 0;
    State magicalState = getModel().magicalState();
    while (states[l] != magicalState) {
      ++l;
    }
      
    return cell.scores[l];
  }
  
  public double backward(SymbolList[] seqs, ScoreType scoreType) 
  throws IllegalSymbolException,
  IllegalAlphabetException,
  IllegalTransitionException {
    return backwardMatrix(seqs, scoreType).getScore();
  }

  public DPMatrix backwardMatrix(SymbolList[] seqs, ScoreType scoreType) 
  throws IllegalSymbolException,
  IllegalAlphabetException,
  IllegalTransitionException {
    if (seqs.length != 2) {
      throw new IllegalArgumentException("This DP object only runs on pairs.");
    }
    lockModel();
    PairDPMatrix matrix = new PairDPMatrix(this, seqs[0], seqs[1]);
    PairDPCursor cursor = new BackMatrixPairDPCursor(
      seqs[0], seqs[1],
      2, 2,
      matrix,
      getEmissionCache(scoreType)
    );
    CellCalculator cc = ccFactory.backwards(scoreType);
    double score = runFB(cursor, cc);
    unlockModel();
    matrix.setScore(score);
    return matrix;
  }

  public DPMatrix backwardMatrix(SymbolList[] seqs, DPMatrix d, ScoreType scoreType) 
  throws IllegalSymbolException,
  IllegalAlphabetException,
  IllegalTransitionException {
    return backwardMatrix(seqs, scoreType);
  }

  public double forward(SymbolList[] seqs, ScoreType scoreType) 
  throws IllegalSymbolException,
  IllegalAlphabetException,
  IllegalTransitionException {
    if (seqs.length != 2) {
      throw new IllegalArgumentException("This DP object only runs on pairs.");
    }
    lockModel();
    PairDPCursor cursor = new LightPairDPCursor(
      seqs[0], seqs[1],
      2, 2, getStates().length, getEmissionCache(scoreType)
    );
    CellCalculator cc = ccFactory.forwards(scoreType);
    double score = runFB(cursor, cc);
    unlockModel();
    return score;
  }

  public DPMatrix forwardMatrix(SymbolList[] seqs, ScoreType scoreType) 
  throws
    IllegalSymbolException,
    IllegalAlphabetException,
    IllegalTransitionException
  {
    if (seqs.length != 2) {
      throw new IllegalArgumentException("This DP object only runs on pairs.");
    }
    lockModel();
    PairDPMatrix matrix = new PairDPMatrix(this, seqs[0], seqs[1]);
    PairDPCursor cursor = new MatrixPairDPCursor(
      seqs[0], seqs[1],
      2, 2, matrix, getEmissionCache(scoreType)
    );
    CellCalculator cc = ccFactory.forwards(scoreType);
    double score = runFB(cursor, cc);
    matrix.setScore(score);
    unlockModel();
    return matrix;
  }

  public DPMatrix forwardMatrix(SymbolList[] seqs, DPMatrix d, ScoreType scoreType) 
  throws
    IllegalSymbolException,
    IllegalAlphabetException,
    IllegalTransitionException
  {
    return forwardMatrix(seqs, scoreType);
  }

  public StatePath viterbi(SymbolList[] seqs, ScoreType scoreType) 
  throws
    IllegalSymbolException,
    IllegalAlphabetException,
    IllegalTransitionException
  {
    if (seqs.length != 2) {
      throw new IllegalArgumentException("This DP object only runs on pairs.");
    }
    lockModel();
    SymbolList seq0 = seqs[0];
    SymbolList seq1 = seqs[1];
    State magic = getModel().magicalState();
    BackPointer TERMINAL_BP = new BackPointer(magic);
    PairDPCursor cursor = new LightPairDPCursor(
      seq0, seq1,
      2, 2, getStates().length, getEmissionCache(scoreType)
    );
    CellCalculator cc = ccFactory.viterbi(scoreType, TERMINAL_BP);
    Cell currentCell = run(cursor, cc);
  
    // Terminate!

    int l = 0;
    State [] states = getStates();
    State magicalState = getModel().magicalState();
    while (states[l] != magicalState) {
      ++l;
    }

    // Traceback...  

    BackPointer[] bpCol = currentCell.backPointers;
    BackPointer bp = bpCol[l];
    List statel = new ArrayList();
    GappedSymbolList gap0 = new SimpleGappedSymbolList(seq0);
    GappedSymbolList gap1 = new SimpleGappedSymbolList(seq1);
    int i0 = seq0.length()+1;
    int i1 = seq1.length()+1;
  
    // parse 1
    //System.out.println("Parse 1");
    for(BackPointer bpi = bp.back; bpi != TERMINAL_BP; bpi = bpi.back) {
      try {
        //System.out.println("bp.back" + bp.back);
      /*System.out.print(
        "Processing " + bpi.state.getName()
      );*/
      statel.add(bpi.state);
      if(bpi.state instanceof EmissionState) { 
        int [] advance = ((EmissionState) bpi.state).getAdvance();
        //System.out.print( "\t" + advance[0] + " " + advance[1]);
        if(advance[0] == 0) {
          gap0.addGapInSource(i0);
          //System.out.println(gap0.seqString());
          //System.out.print("\t-");
        } else {
          i0--;
    	    //System.out.print("\t" + seq0.symbolAt(i0).getToken());
        }
        if(advance[1] == 0) {
          gap1.addGapInSource(i1);
          //System.out.println(gap1.seqString());
          //System.out.print(" -");
        } else {
          i1--;
    	    //System.out.print(" " + seq1.symbolAt(i1).getToken());
        }
      }
      //System.out.println("\tat " + i0 + ", " + i1);
      } catch (IndexOutOfBoundsException ie) {
        while(bpi != TERMINAL_BP) {
          //System.out.println(bpi.state.getName());
          bpi = bpi.back;
        }
        throw new BioError(ie); 
      }
    }
    //System.out.println(gap0.seqString());
    //System.out.println(gap1.seqString());
    double [] scoreA = new double[statel.size()];
    Map aMap = new SmallMap();
    aMap.put(seq0, gap0);
    aMap.put(seq1, gap1);
    Alignment ali = new SimpleAlignment(aMap);
    GappedSymbolList gappedAli = new SimpleGappedSymbolList(ali);

    // parse 2
    //System.out.println("Parse 2");
    int di = statel.size()-1;
    int dj = ali.length()+1;
    for(BackPointer bpi = bp.back; bpi != TERMINAL_BP; bpi = bpi.back) {
      scoreA[di] = bpi.score;
      if(bpi.state instanceof EmissionState) {
        dj--;
      } else {
        gappedAli.addGapInSource(dj);
      }
      di--;
    }

    Collections.reverse(statel);
    SymbolList statesSL = new SimpleSymbolList(getModel().stateAlphabet(), statel);
    SymbolList scoresSL = DoubleAlphabet.fromArray(scoreA);
    StatePath sp = new SimpleStatePath(currentCell.scores[l], gappedAli, statesSL, scoresSL);
    unlockModel();
    return sp;
  }
}
