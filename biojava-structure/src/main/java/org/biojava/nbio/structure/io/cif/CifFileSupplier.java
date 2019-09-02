package org.biojava.nbio.structure.io.cif;

import org.rcsb.cif.model.CifFile;

/**
 * Create a CifFile instance for a given container of structure data.
 * @param <S> the container type used as source
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 * @since 5.3.0
 */
interface CifFileSupplier<S> {
    /**
     * Convert some model instance describing structure information to a CifFile instance.
     * @param container the source of structure information
     * @return a flat CifFile instance, ready for IO operations
     */
    CifFile get(S container);
}
