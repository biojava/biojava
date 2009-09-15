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


package org.biojava.bio.dp;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioError;
import org.biojava.bio.dist.DistributionFactory;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.utils.ChangeVetoException;

/**
 * @author Matthew Pocock
 * @author Lachlan Coin
 */
public class ProfileHMM extends SimpleMarkovModel {
  /**
   * The advance array.
   */
  private final static int [] advance = { 1 };

  /**
   * The number of columns in this model.
   */
  private final int columns;

  /**
   * Match states array.
   * <p>
   * matchStates[0] == matchStates[columns+1] == magicalState().
   */
  private final EmissionState [] matchStates;

  /**
   * Insert states array.
   * <p>
   * From 0 .. columns().
   */
  private final EmissionState [] insertStates;

  /**
   * Delete states array.
   * <p>
   * From 0 .. columns()-1 corresponding to indexes 1..columns().
   */
  private final DotState [] deleteStates;

  /**
   * Retrieve the number of columns in the model.
   *
   * @return the number of columns
   */
  public int columns() {
    return columns;
  }

  /**
   * Retrieve the match state at column indx.
   * <p>
   * The first match state is at index 1, and the last match state is at
   * column columns(). The states at index 0 and columns()+1 are both the
   * magical state. This is so that the whole model backbone can be addressed
   * without writing lots of special-case code.
   *
   * @param indx  the index of the column to retrieve the match state for
   * @return the match state for column indx
   * @throws IndexOutOfBoundsException if indx is negative or above columns()+1
   */
  public EmissionState getMatch(int indx)
  throws IndexOutOfBoundsException {
    if(indx < 0 || indx > (columns+1) ) {
      throw new IndexOutOfBoundsException(
        "Match-state index must be within (0.." + (columns + 1) + "), not " + indx
      );
    }

    return matchStates[indx];
  }

  /**
   * Retrieves the insert state at column indx.
   * <p>
   * Insert_0 is the insert that is accessible directly from the magical state.
   * Insert_1..columns() are 'above' each propper match state. There is no
   * insert state above the magical state at the end of the model, as insert_columns
   * already models trailing inserts.
   *
   * @param indx  the index of the column to retrieve the insert state for
   * @return the insert state for column indx
   * @throws IndexOutOfBoundsException if indx is negative or above columns()
   */
  public EmissionState getInsert(int indx)
  throws IndexOutOfBoundsException {
    if(indx < 0 || indx > columns ) {
      throw new IndexOutOfBoundsException(
        "Insert-state index must be within (0.." + columns + "), not " + indx
      );
    }

    return insertStates[indx];
  }

  /**
   * Retrieves the delete state for column indx.
   * <p>
   * Delete states are 'above' the match state that they delete. There is no
   * delete state for the magical state at either the beginning or end of the
   * model, so the delete state indx can range within (1..columns()).
   *
   * @param indx  the index of the column to retrieve the insert state for
   * @return the insert state for column indx
   * @throws IndexOutOfBoundsException if indx is negative or above columns()
   */
  public DotState getDelete(int indx)
  throws IndexOutOfBoundsException {
    if(indx < 1 || indx > columns ) {
      throw new IndexOutOfBoundsException(
        "delete-state index must be within (1.." + columns + "), not " + indx
      );
    }

    return deleteStates[indx-1];
  }

    /**
     * @deprecated
     */
  public ProfileHMM(
    Alphabet alpha, int columns,
    DistributionFactory matchFactory, DistributionFactory insertFactory) throws IllegalSymbolException, IllegalTransitionException,  IllegalAlphabetException {
    this(alpha, columns, matchFactory, insertFactory,"");
  }



  /**
   * Create a new ProfileHMM.
   * <p>
   * The profile will be over the Alphabet alpha. It will have 'columns' match states,
   * 'columns' delete states and 'columns'+1 insert states. The match states will be
   * created by matchFactory, and the insert states will be created by
   * insertFactory. This gives you great freedom for changing how the states are
   * implemented. For example, insertFactory may always return views onto a single
   * underlying probability distribution in which case all insert states will
   * emit the same type of stuff - or it may create an individual state for
   * each insert state, in which case each insert can be a different type of thing.
   * You could make the insert state a view onto your null model, or onto a
   * protein-specific distribution.
   *
   * @param alpha the Alphabet that the profile hmm will emit
   * @param columns the number of match states
   * @param matchFactory the StateFactory to use for creating the match states
   * @param insertFactory the stateFactory to use for creating the insert states
   */
  public ProfileHMM(
    Alphabet alpha, int columns,
    DistributionFactory matchFactory, DistributionFactory insertFactory, String name)
      throws IllegalSymbolException, IllegalTransitionException,
  IllegalAlphabetException {
    super(1, alpha,name);

    try {
      this.columns = columns;
      this.matchStates = new EmissionState[columns+2];
      this.insertStates = new EmissionState[columns+1];
      this.deleteStates = new DotState[columns];

      EmissionState mO = magicalState();
      EmissionState iO = new SimpleEmissionState(
        "i-0",
        Annotation.EMPTY_ANNOTATION,
        advance,
        insertFactory.createDistribution(alpha)
      );

      matchStates[0] = mO;
      insertStates[0] = iO;

      // first column - a leading insert & the magical state
      addState(iO);

      // 'body' columns
      for(int i = 1; i <= columns; i++) {
        EmissionState mN = new SimpleEmissionState(
          "m-" + i,
          Annotation.EMPTY_ANNOTATION,
          advance,
          matchFactory.createDistribution(alpha)
        );
        EmissionState iN = new SimpleEmissionState(
          "i-" + i,
          Annotation.EMPTY_ANNOTATION,
          advance,
          insertFactory.createDistribution(alpha)
        );
        DotState dN = new SimpleDotState("d-" + i);

        addState(mN);
        addState(iN);
        addState(dN);

        matchStates[i] = mN;
        insertStates[i] = iN;
        deleteStates[i-1] = dN;


        mO = mN;
        iO = iN;
      }
      matchStates[columns+1] = magicalState();
      connectModel();
    } catch (ChangeVetoException cve) {
      throw new BioError(

        "Unable to construct profile HMM", cve
      );
    }
  }

    /** This is called by constructor in setting up the allowed transitions in the model */
    protected void connectModel() throws
        ChangeVetoException, IllegalSymbolException, IllegalTransitionException,IllegalAlphabetException{
        EmissionState mO = getMatch(0);
        EmissionState iO = getInsert(0);
        DotState dO = null;
        createTransition(mO,iO);
        createTransition(iO,iO);
        for(int i = 1; i <= columns(); i++){
            EmissionState mN = getMatch(i);
            EmissionState iN = getInsert(i);
            DotState dN = getDelete(i);

            // from a model state
            createTransition(mO, mN);
            createTransition(mN, iN);
            createTransition(mO, dN);

            // from an insert state
            createTransition(iN, iN);
            createTransition(iO, mN);
            createTransition(iO, dN);

            // from a delete state
            if(i > 1) {
                createTransition(dO, dN);
                createTransition(dO, mN);
            }
            createTransition(dN,iN);

            mO = mN;
            iO = iN;
            dO = dN;
        }
        // for the transitions to end
        createTransition(mO, magicalState());
        createTransition(iO, magicalState());
        createTransition(dO, magicalState());
    }

}
