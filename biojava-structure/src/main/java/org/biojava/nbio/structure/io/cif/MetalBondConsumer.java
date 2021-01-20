package org.biojava.nbio.structure.io.cif;

import org.biojava.nbio.structure.chem.MetalBondDistance;
import org.rcsb.cif.model.Category;

import java.util.List;
import java.util.Map;

public interface MetalBondConsumer extends CifFileConsumer<Map<String,List<MetalBondDistance>>> {
    void consume(Category category);
}
