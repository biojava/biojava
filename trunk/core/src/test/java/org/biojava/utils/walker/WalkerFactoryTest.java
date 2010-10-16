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

package org.biojava.utils.walker;

import junit.framework.TestCase;

import org.biojava.bio.BioException;

/**
 * Tests just the basics of loading classes for viewers. Doesn't use the
 * resulting walkers in any way.
 *
 * @author Matthew Pocock
 */
public class WalkerFactoryTest
extends TestCase {
  public void testNullVisitor() {
    try {
      WalkerFactory.getInstance().getWalker(new Visitor() {});
    } catch (BioException be) {
      throw (AssertionError) new AssertionError(
              "Could not instantiate null visitor").initCause(be);
    }
  }

  public void testAllVisitorNoReturn() {
    try {
      WalkerFactory.getInstance().getWalker(new Visitor() {});
    } catch (BioException be) {
      throw (AssertionError) new AssertionError(
              "Could not instantiate all visitor with no return").initCause(be);
    }
  }


  public void testAllVisitorWithReturn() {
    try {
      WalkerFactory.getInstance().getWalker(new Visitor() {
      });
    } catch (BioException be) {
      throw (AssertionError) new AssertionError(
              "Could not instantiate all visitor with return").initCause(be);
    }
  }

  public void testRepetition() {
    try {
      Visitor vis1 = new Visitor() {};
      Visitor vis2 = new Visitor() {};
      Walker w1 = WalkerFactory.getInstance().getWalker(vis1);
      Walker w2 = WalkerFactory.getInstance().getWalker(vis1);
      Walker w3 = WalkerFactory.getInstance().getWalker(vis2);

      assertEquals("Same class for same visitor: ", w1.getClass(), w2.getClass());
      assertNotSame("Different class for differnt visitor: ", w1.getClass(), w3.getClass());
    } catch (BioException be) {
      throw (AssertionError) new AssertionError(
              "Could not instantiate null visitor").initCause(be);
    }
  }
}
