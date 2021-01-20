package org.biojava.nbio.structure.io.cif;

import org.biojava.nbio.structure.chem.ChemCompDescriptor;
import org.biojava.nbio.structure.chem.ChemicalComponentDictionary;
import org.rcsb.cif.schema.mm.ChemComp;
import org.rcsb.cif.schema.mm.ChemCompAtom;
import org.rcsb.cif.schema.mm.ChemCompBond;

public interface ChemCompConsumer extends CifFileConsumer<ChemicalComponentDictionary> {
    void consumeChemComp(ChemComp c);

    void consumeChemCompAtom(ChemCompAtom atom);

    void consumeChemCompBond(ChemCompBond bond);
}

