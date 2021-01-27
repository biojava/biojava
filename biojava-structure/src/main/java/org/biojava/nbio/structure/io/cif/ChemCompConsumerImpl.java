package org.biojava.nbio.structure.io.cif;

import org.biojava.nbio.structure.chem.ChemicalComponentDictionary;
import org.rcsb.cif.schema.mm.ChemComp;
import org.rcsb.cif.schema.mm.ChemCompAtom;
import org.rcsb.cif.schema.mm.ChemCompBond;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Consumes a CCD file to create the {@link ChemicalComponentDictionary}.
 * @author Sebastian Bittrich
 */
public class ChemCompConsumerImpl implements ChemCompConsumer {
    private static final Logger logger = LoggerFactory.getLogger(ChemCompConsumerImpl.class);
    private final ChemicalComponentDictionary dictionary;
    private String latestChemCompId;

    public ChemCompConsumerImpl() {
        this.dictionary = new ChemicalComponentDictionary();
    }

    @Override
    public void consumeChemComp(ChemComp c) {
        org.biojava.nbio.structure.chem.ChemComp chemComp = new org.biojava.nbio.structure.chem.ChemComp();
        chemComp.setId(c.getId().get(0));
        chemComp.setName(c.getName().get(0));
        chemComp.setType(c.getType().get(0));
        chemComp.setPdbxType(c.getPdbxType().get(0));
        chemComp.setFormula(c.getFormula().get(0));
        chemComp.setMonNstdParentCompId(c.getMonNstdParentCompId().get(0));
        chemComp.setPdbxSynonyms(c.getPdbxSynonyms().get(0));
        chemComp.setPdbxFormalCharge(c.getPdbxFormalCharge().get(0));
        chemComp.setPdbxInitialDate(c.getPdbxInitialDate().get(0));
        chemComp.setPdbxModifiedDate(c.getPdbxModifiedDate().get(0));
        chemComp.setPdbxAmbiguousFlag(c.getPdbxAmbiguousFlag().get(0));
        chemComp.setPdbxReleaseStatus(c.getPdbxReleaseStatus().get(0));
        chemComp.setPdbxReplacedBy(c.getPdbxReplacedBy().get(0));
        chemComp.setPdbxReplaces(c.getPdbxReplaces().get(0));
        chemComp.setFormulaWeight(c.getFormulaWeight().get(0));
        chemComp.setOneLetterCode(c.getOneLetterCode().get(0));
        chemComp.setThreeLetterCode(c.getThreeLetterCode().get(0));
        chemComp.setPdbxModelCoordinatesDetails(c.getPdbxModelCoordinatesDetails().get(0));
        chemComp.setPdbxModelCoordinatesMissingFlag(c.getPdbxModelCoordinatesMissingFlag().get(0));
        chemComp.setPdbxIdealCoordinatesDetails(c.getPdbxIdealCoordinatesDetails().get(0));
        chemComp.setPdbxIdealCoordinatesMissingFlag(c.getPdbxIdealCoordinatesMissingFlag().get(0));
        chemComp.setPdbxModelCoordinatesDbCode(c.getPdbxModelCoordinatesDbCode().get(0));
        chemComp.setPdbxSubcomponentList(c.getPdbxSubcomponentList().get(0));
        chemComp.setPdbxProcessingSite(c.getPdbxProcessingSite().get(0));
        if (chemComp.getId() == null) {
            logger.warn("chem comp ID == null {}", c);
        }
        latestChemCompId = chemComp.getId();
        dictionary.addChemComp(chemComp);
    }

    @Override
    public void consumeChemCompAtom(ChemCompAtom atom) {
        for (int i = 0; i < atom.getRowCount(); i++) {
            org.biojava.nbio.structure.chem.ChemCompAtom a = new org.biojava.nbio.structure.chem.ChemCompAtom();
            a.setCompId(atom.getCompId().get(i));
            a.setAtomId(atom.getAtomId().get(i));
            a.setAltAtomId(atom.getAltAtomId().get(i));
            a.setTypeSymbol(atom.getTypeSymbol().get(i));
            a.setCharge(atom.getCharge().get(i));
            a.setPdbxAlign(atom.getPdbxAlign().get(i));
            a.setPdbxAromaticFlag(atom.getPdbxAromaticFlag().get(i));
            a.setPdbxLeavingAtomFlag(atom.getPdbxLeavingAtomFlag().get(i));
            a.setPdbxStereoConfig(atom.getPdbxStereoConfig().get(i));
            a.setModelCartnX(atom.getModelCartnX().get(i));
            a.setModelCartnY(atom.getModelCartnY().get(i));
            a.setModelCartnZ(atom.getModelCartnZ().get(i));
            a.setPdbxModelCartnXIdeal(atom.getPdbxModelCartnXIdeal().get(i));
            a.setPdbxModelCartnYIdeal(atom.getPdbxModelCartnYIdeal().get(i));
            a.setPdbxModelCartnZIdeal(atom.getPdbxModelCartnZIdeal().get(i));
            a.setPdbxComponentAtomId(atom.getPdbxComponentAtomId().get(i));
            a.setPdbxComponentCompId(atom.getPdbxComponentCompId().get(i));
            a.setPdbxOrdinal(atom.getPdbxOrdinal().get(i));
            dictionary.getChemComp(latestChemCompId).getAtoms().add(a);
        }
    }

    @Override
    public void consumeChemCompBond(ChemCompBond bond) {
        for (int i = 0; i < bond.getRowCount(); i++) {
            org.biojava.nbio.structure.chem.ChemCompBond b = new org.biojava.nbio.structure.chem.ChemCompBond();
            b.setAtomId1(bond.getAtomId1().get(i));
            b.setAtomId2(bond.getAtomId2().get(i));
            b.setCompId(bond.getCompId().get(i));
            b.setPdbxAromaticFlag(bond.getPdbxAromaticFlag().get(i));
            b.setPdbxOrdinal(bond.getPdbxOrdinal().get(i));
            b.setPdbxStereoConfig(bond.getPdbxStereoConfig().get(i));
            b.setValueOrder(bond.getValueOrder().get(i));
            dictionary.getChemComp(latestChemCompId).getBonds().add(b);
        }
    }

    @Override
    public void prepare() {

    }

    @Override
    public void finish() {

    }

    @Override
    public ChemicalComponentDictionary getContainer() {
        return dictionary;
    }
}
