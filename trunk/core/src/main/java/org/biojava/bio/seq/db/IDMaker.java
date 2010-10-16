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

package org.biojava.bio.seq.db;

import java.io.IOException;
import java.io.Serializable;

import org.biojava.bio.BioError;
import org.biojava.bio.seq.Sequence;
import org.biojava.utils.StaticMemberPlaceHolder;

/**
 * Interface for objects that define how to make an ID for a sequence.
 * <p>
 * Nine times out of ten, you will use one of HashSequenceDB.byURN or
 * HashSequenceDB.byName, but once in a blue-moon, you will want some other
 * systematic way of retrieveing Sequences. This interface is here to allow
 * you to plug in this functionality if you need it.
 *
 * @author Matthew Pocock
 */
public interface IDMaker {
  /**
   * Calculate the id for a sequence.
   * <p>
   * Each unique sequence should return a unique ID.
   *
   * @param seq the sequence to ID
   * @return the id for the sequence
   */
  String calcID(Sequence seq);


  /**
   * A simple implementation of IDMaker that hashes by URN.
   *
   */
  public final static IDMaker byURN = new ByURN();

  static class ByURN implements Serializable, IDMaker {
    public String calcID(Sequence seq) {
      return seq.getURN();
    }
    private Object writeReplace() throws IOException {
      try {
        return new StaticMemberPlaceHolder(
          IDMaker.class.getField("byURN")
        );
      } catch (NoSuchFieldException nsfe) {
        throw new BioError(
          "Could not find field while serializing",
          nsfe
        );
      }
    }
  }

  /**
   * A simple implementation of IDMaker that hashes by sequence name.
   *
   */
  public final static IDMaker byName = new ByName();

  static class ByName implements Serializable, IDMaker {
    public String calcID(Sequence seq) {
      return seq.getName();
    }
    private Object writeReplace() throws IOException {
      try {
        return new StaticMemberPlaceHolder(
          IDMaker.class.getField("byName")
        );
      } catch (NoSuchFieldException nsfe) {
        throw new BioError(
          "Could not find field while serializing",
          nsfe
        );
      }
    }
  }
}

