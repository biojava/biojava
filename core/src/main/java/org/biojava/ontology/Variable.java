package org.biojava.ontology;

/**
 *
 *
 * @author Matthew Pocock
 */
public interface Variable
extends Term {
  public static class Impl
          extends Term.Impl
          implements Variable
 {
    public Impl(Ontology ontology, String name, String description) {
      super(ontology, name, description);
    }
    public Impl(Ontology ontology, String name, String description, Object[] synonyms) {
      super(ontology, name, description, synonyms);
    }
  }
}
