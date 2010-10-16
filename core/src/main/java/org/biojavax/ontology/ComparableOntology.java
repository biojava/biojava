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

package org.biojavax.ontology;

import java.util.Set;

import org.biojava.ontology.Ontology;
import org.biojava.ontology.Term;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.Changeable;

/**
 * An Ontology that can be compared to another.
 * @author Richard Holland
 * @see ComparableTerm
 * @see ComparableTriple
 * @since 1.5
 */
public interface ComparableOntology extends Ontology,Comparable,Changeable {
    
    public static final ChangeType TERM = new ChangeType(
            "This ontology's terms have changed",
            "org.biojavax.ontology.ComparableOntology",
            "TERM"
            );
    public static final ChangeType TRIPLE = new ChangeType(
            "This ontology's triples have changed",
            "org.biojavax.ontology.ComparableOntology",
            "TRIPLE"
            );
    public static final ChangeType DESCRIPTION = new ChangeType(
            "This ontology's description has changed",
            "org.biojavax.ontology.ComparableOntology",
            "DESCRIPTION"
            );
    
    /**
     * Sets a human-readable description of this ontology.
     * @param description the description.
     * @throws ChangeVetoException in case of problems.
     */
    public void setDescription(String description) throws ChangeVetoException;
    
    /**
     * Return a human-readable description of this ontology.
     * @return the description.
     */
    public String getDescription();
    
    /**
     * Clears out all the terms and populates the ontology with the contents
     * of the set passed. The terms should be ComparableTerms.
     * @param terms a set of Term objects this ontology should have.
     * @throws ChangeVetoException if any of them are unacceptable.
     * @see ComparableTerm
     */
    public void setTermSet(Set terms) throws ChangeVetoException;
    
    /**
     * Returns the set of terms in this ontology.
     * @return a set of ComparableTerm objects.
     * @see ComparableTerm
     */
    public Set getTermSet();
    
    /**
     * Clears out all the triples and populates the ontology with the contents
     * of the set passed.
     * @param triples the set of ComparableTriple objects this ontology should have.
     * @throws ChangeVetoException if any of them are unacceptable.
     * @see ComparableTriple
     */
    public void setTripleSet(Set triples) throws ChangeVetoException;
    
    /**
     * Returns the set of triples in this ontology.
     * @return the set of ComparableTriple objects.
     */
    public Set getTripleSet();
    
    /**
     * Looks for a term with the given name and returns it. If it couldn't be found,
     * then it creates it, adds it to the ontology, then returns it.
     * @param name the name of the term to look for.
     * @return the ComparableTerm representing that name.
     */
    public ComparableTerm getOrCreateTerm(String name);
    
    /**
     * Looks for a triple with the given subject object and predicate and returns it. 
     * If it couldn't be found, then it creates it, adds it to the ontology,
     * then returns it.
     * @param subject the subject of the triple eg apple
     * @param object the object of the triple eg fruit
     * @param predicate the relationship of the triple eg is_a
     * @return the ComparableTriple representing the object subject and predicate.
     */
    public ComparableTriple getOrCreateTriple(Term subject, Term object, Term predicate);
    /**
     * Looks for a term with the same name as the given term and returns it.
     * If it couldn't be found, then it creates it, adds it to the ontology,
     * then returns it.
     * @param term the term to look for.
     * @return the ComparableTerm representing that term in this ontology.
     */
    public ComparableTerm getOrImportTerm(Term term);
}

