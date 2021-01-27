package org.biojava.nbio.structure.io.cif;

import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.rcsb.cif.model.CifFile;

import java.util.ArrayList;
import java.util.List;

/**
 * Convert a structure to a CifFile.
 * @author Sebastian Bittrich
 */
public class CifStructureSupplierImpl extends AbstractCifFileSupplier<Structure> {
    @Override
    public CifFile get(Structure container) {
        return getInternal(container, collectWrappedAtoms(container));
    }

    private List<WrappedAtom> collectWrappedAtoms(Structure structure) {
        List<WrappedAtom> wrappedAtoms = new ArrayList<>();

        for (int modelIndex = 0; modelIndex < structure.nrModels(); modelIndex++) {
            final int model = modelIndex + 1;
            for (Chain chain : structure.getChains(modelIndex)) {
                handleChain(chain, model, wrappedAtoms);
            }
        }

        return wrappedAtoms;
    }
}
