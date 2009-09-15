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
 */

package org.biojava.bio.seq.db;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.NoSuchElementException;
import java.util.Set;

import junit.framework.TestCase;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.SequenceTools;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.ChangeVetoException;


/**
 * Tests for HashSequenceDB
 * @author Mark Schreiber
 */
public class TestHashSequenceDB extends TestCase {
  private HashSequenceDB hashSequenceDB = null;
  private HashSequenceDB hashSequenceDB2 = null;
  private Sequence seq = null;

  public TestHashSequenceDB(String name) {
    super(name);
  }

  protected void setUp() throws Exception {
    super.setUp();
    hashSequenceDB = new HashSequenceDB();
    hashSequenceDB2 = new HashSequenceDB(new IDMaker.ByURN(), "seqDB");
    SymbolList syms = ProteinTools.createProtein("hgfds");
    Sequence seq = SequenceTools.createSequence(syms, "urn:biojava:bar", "bar",
                                              Annotation.EMPTY_ANNOTATION);

    hashSequenceDB.addSequence(seq);
    hashSequenceDB2.addSequence(seq);

  }

  protected void tearDown() throws Exception {
    hashSequenceDB = null;
    super.tearDown();
  }



  public void testAddSequence1() throws ChangeVetoException {
    String id = "myseq";
    hashSequenceDB.addSequence(id, seq);
    hashSequenceDB2.addSequence(id, seq);
  }

  public void testGetName() {
    String expectedReturn = "seqDB";
    String actualReturn = hashSequenceDB2.getName();
    assertEquals("return value", expectedReturn, actualReturn);
  }

  public void testGetSequence() throws IllegalIDException {
    String id = "bar";
    String urn = "urn:biojava:bar";
    Sequence s = hashSequenceDB.getSequence(id);
    if(s == null)
      fail("no sequence returned for :"+id);
    assertTrue(s.getName().equals (id));

    s = hashSequenceDB2.getSequence(urn);
    if(s == null)
      fail("no sequence returned for :"+urn);
    assertTrue(s.getURN().equals(urn));

  }

  public void testIds() {
    Set ids = hashSequenceDB.ids();
    assertNotNull(ids);
    assertTrue(ids.contains("bar"));
    assertTrue(ids.size() == 1);

    ids = hashSequenceDB2.ids();
    assertNotNull(ids);
    assertTrue(ids.contains("urn:biojava:bar"));
    assertTrue(ids.size() == 1);
  }

  public void testRemoveSequence() throws BioException, ChangeVetoException {
    String id = "bar";
    String urn = "urn:biojava:bar";
    hashSequenceDB.removeSequence(id);
    hashSequenceDB2.removeSequence(urn);

    assertFalse(hashSequenceDB.ids().contains(id));
    assertFalse(hashSequenceDB.ids().contains(urn));

    assertTrue(hashSequenceDB.ids().size() == 0);
    assertTrue(hashSequenceDB2.ids().size() == 0);
  }

  public void testSequenceIterator() {
    SequenceIterator iter = hashSequenceDB.sequenceIterator();
    assertNotNull(iter);

    Sequence s = null;
    try {
      s = iter.nextSequence();
    } catch (BioException ex) {
      fail(ex.getMessage());
    } catch (NoSuchElementException ex) {
      ex.getMessage();
    }
    assertNotNull(s);
    //should be no more
    assertFalse(iter.hasNext());
  }

  public void testSerialization() throws IOException, ClassNotFoundException{
    ByteArrayOutputStream baos = new ByteArrayOutputStream();
    ObjectOutputStream oos = new ObjectOutputStream(baos);
    oos.writeObject(hashSequenceDB);
    oos.writeObject(hashSequenceDB2);
    oos.flush();

    ByteArrayInputStream bais = new ByteArrayInputStream(baos.toByteArray());
    ObjectInputStream ois = new ObjectInputStream(bais);
    HashSequenceDB seqdb = (HashSequenceDB)ois.readObject();
    HashSequenceDB seqdb2 = (HashSequenceDB)ois.readObject();
    bais.close();
    ois.close();
    oos.close();
    baos.close();

    Sequence s = null;
    try {
      assertNotNull(seqdb);
      s = seqdb.getSequence("bar");
      assertNotNull(s);
    } catch (IllegalIDException ex) {
      fail(ex.getMessage());
    }

    try {
      assertNotNull(seqdb2);
      s = seqdb2.getSequence("urn:biojava:bar");
      assertNotNull(s);
    } catch (IllegalIDException ex1) {
      fail(ex1.getMessage());
    }
  }
}
