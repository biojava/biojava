package org.biojava.nbio.structure.io.cif;

import org.biojava.nbio.structure.chem.MetalBondDistance;
import org.rcsb.cif.model.Category;
import org.rcsb.cif.model.FloatColumn;
import org.rcsb.cif.model.StrColumn;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by andreas on 6/9/16.
 */
public class MetalBondConsumerImpl implements MetalBondConsumer {
    private final Map<String, List<MetalBondDistance>> definitions = new HashMap<>();

    @Override
    public void prepare() {
        definitions.clear();
    }

    @Override
    public void finish() {
        // minimize memory consumption
        for (List<MetalBondDistance> d : definitions.values()){
            ((ArrayList<MetalBondDistance>) d).trimToSize();
        }
    }

    @Override
    public void consume(Category category) {
        StrColumn atomType1 = (StrColumn) category.getColumn("atom_type_1");
        StrColumn atomType2 = (StrColumn) category.getColumn("atom_type_2");
        FloatColumn lowerLimit = (FloatColumn) category.getColumn("lower_limit");
        FloatColumn upperLimit = (FloatColumn) category.getColumn("upper_limit");
        for (int i = 0; i < category.getRowCount(); i++) {
            MetalBondDistance d = new MetalBondDistance();

            d.setAtomType1(atomType1.get(i));
            d.setAtomType2(atomType2.get(i));
            d.setLowerLimit((float) lowerLimit.get(i));
            d.setUpperLimit((float) upperLimit.get(i));

            List<MetalBondDistance> defs = definitions.computeIfAbsent(d.getAtomType1(), k -> new ArrayList<>());
            defs.add(d);
        }
    }

    @Override
    public Map<String, List<MetalBondDistance>> getContainer(){
        return definitions;
    }
}
