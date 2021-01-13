package org.biojava.nbio.structure.chem;

import org.rcsb.cif.schema.mm.ChemComp;

public class ChemCompContainer {
    private final ChemComp delegate;

    public ChemCompContainer(ChemComp chemComp) {
        this.delegate = chemComp;
    }


}
