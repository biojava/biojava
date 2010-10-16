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

/**
 * <code>TermImplTest</code>
 *
 * @author Moses Hohman
 * @since 1.4
 */

public class TermImplTest extends TestCase {
    private Ontology onto;

    protected void setUp() throws Exception {
        onto = OntoTools.getDefaultFactory().createOntology("test", "");
    }

    //no reason to override Object.equals() since in a given ontology all terms should be unique objects
    public void testUnequalWhenNotSameObject() {
        assertFalse(new Term.Impl(onto, new String("test term"), new String("")).equals(new Term.Impl(onto, new String("test term"), new String(""))));
    }

    public void testToString() {
        assertEquals("name", new Term.Impl(onto, "name", "").toString());
    }

    public void testGetAnnotationIsBlank() {
        assertEquals("{}", new Term.Impl(onto, "", "").getAnnotation().toString());
    }

    public void testCannotConstructWithNullName() {
        try {
            new Term.Impl(onto, null, "");
            fail("Should have thrown a NullPointerException");
        } catch (NullPointerException expected) {
        }
    }

    public void testCannotConstructWithNullOntology() {
        try {
            new Term.Impl(null, "", "");
            fail("Should have thrown a NullPointerException");
        } catch (NullPointerException expected) {
        }
    }

    // Description of terms from now on can change!
//    public void testCannotConstructWithNullDescription() {
//        try {
//            new Term.Impl(onto, "", null);
//            fail("Should have thrown a NullPointerException");
//        } catch (NullPointerException expected) {
//        }
//    }
}