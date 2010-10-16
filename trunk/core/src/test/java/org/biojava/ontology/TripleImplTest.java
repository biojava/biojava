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

package org.biojava.ontology;

import junit.framework.TestCase;

import org.biojava.utils.ChangeVetoException;

/**
 * <code>TripleImplTest</code>
 *
 * @author Moses Hohman
 * @author Matthew Pocock
 * @since 1.4
 */

public class TripleImplTest extends TestCase {
    private Ontology onto;
    private Term subject;
    private Term object;
    private Term relation;

    protected void setUp() throws Exception {
        onto = OntoTools.getDefaultFactory().createOntology("test", "");
        subject = onto.createTerm("subject", "");
        object = onto.createTerm("object", "");
        relation = onto.importTerm(OntoTools.IS_A, null);
    }

    public void testEqualsWhenEqual() throws ChangeVetoException, AlreadyExistsException {
        Triple.Impl triple1 = new Triple.Impl(subject, object, relation);
        Triple.Impl triple2 = new Triple.Impl(subject, object, relation);
        assertEquals(triple1, triple2);
    }

    public void testEqualsWhenSubjectsUnequal() throws ChangeVetoException, AlreadyExistsException {
        Triple.Impl triple1 = new Triple.Impl(subject, object, relation);
        Triple.Impl triple2 = new Triple.Impl(onto.createTerm("something else", ""), object, relation);
        assertFalse(triple1.equals(triple2));
    }

    public void testEqualsWhenObjectsUnequal() throws ChangeVetoException, AlreadyExistsException {
        Triple.Impl triple1 = new Triple.Impl(subject, object, relation);
        Triple.Impl triple2 = new Triple.Impl(subject, onto.createTerm("something else", ""), relation);
        assertFalse(triple1.equals(triple2));
    }

    public void testEqualsWhenRelationsUnequal() throws ChangeVetoException, AlreadyExistsException {
        Triple.Impl triple1 = new Triple.Impl(subject, object, relation);
        Triple.Impl triple2 = new Triple.Impl(subject, object, onto.importTerm(OntoTools.RELATION, null));
        assertFalse(triple1.equals(triple2));
    }

    public void testHashcodesEqualWhenEqual() throws ChangeVetoException, AlreadyExistsException {
        Triple.Impl triple1 = new Triple.Impl(subject, object, relation);
        Triple.Impl triple2 = new Triple.Impl(subject, object, relation);
        assertEquals(triple1.hashCode(), triple2.hashCode());
    }

    public void testCannotConstructWithNullSubject() throws ChangeVetoException, AlreadyExistsException {
        try {
            new Triple.Impl(null, object, relation);
            fail("should have thrown an IllegalArgumentException");
        } catch (NullPointerException expected) {
        }
    }

    public void testCannotConstructWithNullObject() {
        try {
            new Triple.Impl(subject, null, relation);
            fail("should have thrown an IllegalArgumentException");
        } catch (NullPointerException expected) {
        }
    }

    public void testCannotConstructWithNullRelation() {
        try {
            new Triple.Impl(subject, object, null);
            fail("should have thrown an IllegalArgumentException");
        } catch (NullPointerException expected) {
        }
    }
}