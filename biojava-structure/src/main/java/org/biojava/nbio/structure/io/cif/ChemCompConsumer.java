package org.biojava.nbio.structure.io.cif;

import org.biojava.nbio.structure.chem.ChemicalComponentDictionary;
import org.rcsb.cif.schema.mm.ChemComp;
import org.rcsb.cif.schema.mm.ChemCompAtom;
import org.rcsb.cif.schema.mm.ChemCompBond;

/**
 * Create the {@link ChemicalComponentDictionary} from CIF data.
 * @author Sebastian Bittrich
 * @since 6.0.0
 */
public interface ChemCompConsumer extends CifFileConsumer<ChemicalComponentDictionary> {
    /**
     * Consume a particular Cif category.
     * @param c data
     */
    void consumeChemComp(ChemComp c);

    /**
     * Consume a particular Cif category.
     * @param atom data
     */
    void consumeChemCompAtom(ChemCompAtom atom);

    /**
     * Consume a particular Cif category.
     * @param bond data
     */
    void consumeChemCompBond(ChemCompBond bond);
}

