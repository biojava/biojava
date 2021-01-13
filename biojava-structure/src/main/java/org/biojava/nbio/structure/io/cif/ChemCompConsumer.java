package org.biojava.nbio.structure.io.cif;

import org.biojava.nbio.structure.chem.ChemComp;
import org.rcsb.cif.schema.mm.ChemCompAtom;
import org.rcsb.cif.schema.mm.ChemCompBond;
import org.rcsb.cif.schema.mm.PdbxChemCompDescriptor;

public interface ChemCompConsumer extends CifFileConsumer<ChemComp> {
    void consumeChemCompAtom(ChemCompAtom chemCompAtom);

    void consumeChemCompBond(ChemCompBond chemCompBond);

    void consumePdbxChemCompDescriptor(PdbxChemCompDescriptor pdbxChemCompDescriptor);
}
