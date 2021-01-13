package org.biojava.nbio.structure.io.cif;

import org.rcsb.cif.model.Category;

import java.io.Serializable;

/**
 * Flag for BioJava beans that actually resemble categories defined by the mmCIF schema.
 * @param <C> the modeled ciftools-java category
 */
public interface CifBean<C extends Category> extends Serializable {
}
