package org.biojava.nbio.structure.io.cif;

import org.biojava.nbio.structure.chem.ChemComp;
import org.biojava.nbio.structure.chem.ChemicalComponentDictionary;
import org.biojava.nbio.structure.chem.ResidueType;
import org.rcsb.cif.schema.mm.ChemCompAtom;
import org.rcsb.cif.schema.mm.ChemCompBond;
import org.rcsb.cif.schema.mm.PdbxChemCompDescriptor;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class ChemCompConsumerImpl implements ChemCompConsumer<ChemComp> {
    private static final Logger logger = LoggerFactory.getLogger(ChemCompConsumerImpl.class);
    private ChemicalComponentDictionary dictionary;
    private String latestChemCompId;

    public ChemCompConsumerImpl(){
        dictionary = new ChemicalComponentDictionary();
    }

    @Override
    public void prepare() {

    }

    public ChemicalComponentDictionary getDictionary(){
        return dictionary;
    }

    @Override
    public void consumeChemCompAtom(ChemCompAtom chemCompAtom) {
        dictionary.getChemComp(latestChemCompId).getAtoms().add(chemCompAtom);
    }

    @Override
    public void consumeChemCompBond(ChemCompBond chemCompBond) {
        dictionary.getChemComp(latestChemCompId).getBonds().add(chemCompBond);
    }

    @Override
    public void consumePdbxChemCompDescriptor(PdbxChemCompDescriptor pdbxChemCompDescriptor) {
        ChemComp cc = dictionary.getChemComp(latestChemCompId);
        cc.getDescriptors().add(pdbxChemCompDescriptor);
    }

    @Override
    public void finish() {

    }

    @Override
    public ChemComp getContainer() {
        if (c.getId() == null)
            logger.warn("chem comp ID == null " + c);

        latestChemCompId = c.getId();
        dictionary.addChemComp(c);
        if (c.getResidueType() == ResidueType.nonPolymer) {
            return;
        }

        if (c.getResidueType() == ResidueType.saccharide) {
            return;
        }

        if (c.getResidueType() == ResidueType.dSaccharide) {
            return;
        }
    }
}
