package org.biojava.nbio.structure.io.cif;

import org.biojava.nbio.structure.Chain;
import org.rcsb.cif.model.CifFile;

import java.util.ArrayList;
import java.util.List;

/**
 * Convert a chain to a {@link CifFile}.
 * @author Sebastian Bittrich
 */
public class CifChainSupplierImpl extends AbstractCifFileSupplier<Chain> {
    @Override
    public CifFile get(Chain container) {
        return getInternal(container.getStructure(), collectWrappedAtoms(container));
    }

    private List<WrappedAtom> collectWrappedAtoms(Chain chain) {
        List<WrappedAtom> wrappedAtoms = new ArrayList<>();
        handleChain(chain, 1, wrappedAtoms);
        return wrappedAtoms;
    }
}
