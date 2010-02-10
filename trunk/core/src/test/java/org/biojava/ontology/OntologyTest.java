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
 * Tests for Ontology.
 *
 * @author Thomas Down
 * @author Matthew Pocock
 * @author Moses Hohman
 * @since 1.4
 */
public class OntologyTest extends TestCase {
    private String name;
    private String description;
    private Ontology onto;
    private Term animal;
    private Term fish;
    private Term mammal;
    private Term human;
    private Term isa;
    private Triple fish_isa_animal;
    private Triple mammal_isa_animal;
    private Triple human_isa_mammal;

    protected void setUp() throws Exception {
        name = "Tester";
        description = "Our Description";
        onto = OntoTools.getDefaultFactory().createOntology(name, description);

        isa = onto.importTerm(OntoTools.IS_A, null);
        animal = onto.createTerm("Animal", "An animal");
        fish = onto.createTerm("Fish", "A swimming, cold-blooded thingy");
        mammal = onto.createTerm("Mammal", "A milk-producing quadraped");
        human = onto.createTerm("Human", "Us");
        fish_isa_animal = onto.createTriple(fish, animal, isa, null, null);
        mammal_isa_animal = onto.createTriple(mammal, animal, isa, null, null);
        human_isa_mammal = onto.createTriple(human, mammal, isa, null, null);
    }

    // basic properties
    public void testName() {
        assertEquals(onto.getName(), name);
    }

    public void testDescription() {
        assertEquals(onto.getDescription(), description);
    }

    //terms
    public void testTermsSize() {
        assertEquals(8, onto.getTerms().size());
    }

    public void testTermsContainAnimal() {
        assertTrue(onto.getTerms().contains(animal));
    }

    public void testTermsContainFish() {
        assertTrue(onto.getTerms().contains(fish));
    }

    public void testTermsContainMammal() {
        assertTrue(onto.getTerms().contains(mammal));
    }

    public void testTermsContainHuman() {
        assertTrue(onto.getTerms().contains(human));
    }

    //terms by name
    public void testGetAnimalByName() {
        assertEquals(animal, onto.getTerm("Animal"));
    }

    // triples
    public void testGetAllTriplesSize() {
        assertEquals(3, onto.getTriples(null, null, null).size());
    }

    public void testGetAllTriplesContainsFishIsaAnimal() {
        assertTrue(onto.getTriples(null, null, null).contains(fish_isa_animal));
    }

    public void testGetAllTriplesContainsMammalIsaAnimal() {
        assertTrue(onto.getTriples(null, null, null).contains(mammal_isa_animal));
    }

    public void testGetAllTriplesContainsHumanIsaMammal() {
        assertTrue(onto.getTriples(null, null, null).contains(human_isa_mammal));
    }

    // triple searching
    public void testTripleSearchingByObject() {
        assertEquals(2, onto.getTriples(null, animal, null).size());
        assertEquals(1, onto.getTriples(null, mammal, null).size());
        assertEquals(0, onto.getTriples(null, human, null).size());
    }

    public void testTripleSearchBySubject() {
        assertEquals(1, onto.getTriples(onto.getTerm("Human"), null, null).size());
    }

    public void testTripleSearchingByRelation() {
        assertEquals(3, onto.getTriples(null, null, isa).size());
    }

    public void testTripleSearchingWithAllTerms() {
        assertEquals(1, onto.getTriples(human, mammal, isa).size());
        assertEquals(0, onto.getTriples(human, animal, isa).size());
    }

    public void testTripleSearchForSameTermWithDifferentCriteria() {
        assertEquals(onto.getTriples(human, null, isa), onto.getTriples(null, mammal, isa));
    }

    public void testCannotCreateDuplicateTerm() throws ChangeVetoException {
        try {
            onto.createTerm("Human", "Us");
            fail("should have thrown an AlreadyExistsException");
        } catch (AlreadyExistsException expected) {
        }
    }

    public void testCannotCreateTripleTermWithForeignSubject() throws OntologyException, ChangeVetoException {
        try {
            Ontology onto2 = OntoTools.getDefaultFactory().createOntology("other", "");
            onto.createTriple(onto2.createTerm("foreign term", ""), human, isa, null, null);
            fail("Should have thrown an IllegalArgumentException");
        } catch (IllegalArgumentException expected) {
        }
    }

    public void testCannotCreateTripleTermWithForeignObject() throws OntologyException, ChangeVetoException {
        try {
            Ontology onto2 = OntoTools.getDefaultFactory().createOntology("other", "");
            onto.createTriple(fish, onto2.createTerm("foreign term", ""), isa, null, null);
            fail("Should have thrown an IllegalArgumentException");
        } catch (IllegalArgumentException expected) {
        }
    }

    public void testCannotCreateTripleTermWithForeignRelation() throws OntologyException, ChangeVetoException {
        try {
            Ontology onto2 = OntoTools.getDefaultFactory().createOntology("other", "");
            onto.createTriple(human, fish, onto2.createTerm("wants a", ""), null, null);
            fail("Should have thrown an IllegalArgumentException");
        } catch (IllegalArgumentException expected) {
        }
    }
}

