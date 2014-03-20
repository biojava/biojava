package org.biojava3.ontology;

/**
 * A factory for Ontology instances.
 *
 * @author Matthew Pocock
 */

public interface OntologyFactory {
  /**
   * Creates a new Ontology
   *
   * @param name  the name to give the ontology
   * @param description the description for the ontology
   * @return an Ontology 
   * @throws NullPointerException if either name or description are null
   * @throws OntologyException if the ontology could not be created
   */
  public Ontology createOntology(String name, String description)
  throws OntologyException;
}
