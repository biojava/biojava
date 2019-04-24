package org.biojava.nbio.structure.io.cif;

import org.biojava.nbio.structure.Structure;
import org.rcsb.cif.model.CifFile;

/**
 * Convert a BioJava {@link Structure} to a CifFile.
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 * @since 5.2.1
 */
public class CifFileSupplierImpl implements CifFileSupplier<Structure> {
    @Override
    public CifFile get(Structure container) {
        // TODO impl
        return null;
    }
}
