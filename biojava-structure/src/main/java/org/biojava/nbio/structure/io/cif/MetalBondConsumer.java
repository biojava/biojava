package org.biojava.nbio.structure.io.cif;

import org.biojava.nbio.structure.chem.MetalBondDistance;
import org.rcsb.cif.model.Category;

import java.util.List;
import java.util.Map;

/**
 * Consume metal bond data.
 * @author Sebastian Bittrich
 * @since 6.0.0
 */
public interface MetalBondConsumer extends CifFileConsumer<Map<String, List<MetalBondDistance>>> {
    void consume(Category category);
}
