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

import java.io.Serializable;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeVetoException;

/**
 * An implementation of SequenceDB that uses an underlying HashMap to store the
 * sequence objects.
 *
 * @author Matthew Pocock
 * @author <A href="mailto:Gerald.Loeffler@vienna.at">Gerald Loeffler</A>
 */
public class HashSequenceDB
  extends
    AbstractSequenceDB
  implements
    SequenceDB,
    Serializable {
  /**
   * The sequence-by-id map.
   */
  final private Map sequenceByID;
  
  /**
   * An object to extract an ID for a sequence.
   */
  final private org.biojava.bio.seq.db.IDMaker idMaker;

  /** 
   * The name of this sequence database.
   */
  private String name;
  
  public String getName() {
    return name;
  }

  public Sequence getSequence(String id) 
      throws IllegalIDException
  {
      Sequence seq = (Sequence) sequenceByID.get(id);
      if (seq == null) {
          throw new IllegalIDException("Sequence with ID " + id + " could not be found");
      }
      return seq;
  }

  public Set ids() {
    return sequenceByID.keySet();
  }

  public SequenceIterator sequenceIterator() {
    return new SequenceIterator() {
      Iterator seqI = sequenceByID.values().iterator();
      public boolean hasNext() { return seqI.hasNext(); }
      public Sequence nextSequence() { return (Sequence) seqI.next(); }
    };
  }

  /**
   * Add a sequence under a particular id.
   *
   * @param id  the id to use
   * @param seq the Sequence to add
   * @throws ChangeVetoException if this addition was vetoed
   */
  public void addSequence(String id, Sequence seq)
  throws ChangeVetoException {
    if(!hasListeners()) {
      sequenceByID.put(id, seq);
    } else {
      ChangeSupport changeSupport = getChangeSupport(SequenceDB.SEQUENCES);
      synchronized(changeSupport) {
        ChangeEvent ce = new ChangeEvent(
          this,
          SequenceDB.SEQUENCES,
          new Object[] { id, seq },
          null
        );
        changeSupport.firePreChangeEvent(ce);
        sequenceByID.put(id, seq);
        changeSupport.firePostChangeEvent(ce);
      }
    }
  }

  /**
   * Retrieve the IDMaker associated with this database.
   *
   * @return the current IDMaker object
   */
  public org.biojava.bio.seq.db.IDMaker getIDMaker() {
    return idMaker;
  }
  
  public void addSequence(Sequence seq)
  throws ChangeVetoException {
    String id = idMaker.calcID(seq);
    if(!hasListeners()) {
      addSequence(id, seq);
    } else {
      ChangeSupport changeSupport = getChangeSupport(SequenceDB.SEQUENCES);
      synchronized(changeSupport) {
        ChangeEvent ce = new ChangeEvent(
          this,
          SequenceDB.SEQUENCES,
          id,
          null
        );
        changeSupport.firePreChangeEvent(ce);
        sequenceByID.remove(id);
        changeSupport.firePostChangeEvent(ce);
      }
    }
  }
  
  public void removeSequence(String id)
  throws BioException, ChangeVetoException {
    if(!hasListeners()) {
      sequenceByID.remove(id);
    } else {
      ChangeSupport changeSupport = getChangeSupport(SequenceDB.SEQUENCES);
      synchronized(changeSupport) {
        ChangeEvent ce = new ChangeEvent(
          this,
          SequenceDB.SEQUENCES,
          null,
          id
        );
        changeSupport.firePreChangeEvent(ce);
        sequenceByID.remove(id);
        changeSupport.firePostChangeEvent(ce);
      }
    }
  }

  /**
   * Generate a HashSequenceDB object that will use byName to generate ids for
   * sequences and have a null name.
   */
  public HashSequenceDB() {
    this(IDMaker.byName, null);
  }
  
  /**
   * Generate a HashSequenceDB object that will use idMaker to generate ids for
   * sequences and have a null name.
   *
   * @param idMaker the object that will work out the default id for a sequence
   */
  public HashSequenceDB(org.biojava.bio.seq.db.IDMaker idMaker) {
    this(idMaker, null);
  }
  
  /**
   * Generate a HashSequenceDB object that will use byName to generate ids and
   * will have the requested name.
   *
   * @param name the name for this database
   */
  public HashSequenceDB(String name) {
    this(IDMaker.byName, name);
  }
  
  /**
   * Generate a HashSequenceDB object that will use idMaker to generate ids for
   * sequences and have the requested name.
   *
   * @param idMaker the object that will work out the default id for a sequence
   * @param name the name for this database
   */
  public HashSequenceDB(org.biojava.bio.seq.db.IDMaker idMaker, String name) {
    this.idMaker = idMaker;
    this.name = name;
    this.sequenceByID = new HashMap();
  }
}
