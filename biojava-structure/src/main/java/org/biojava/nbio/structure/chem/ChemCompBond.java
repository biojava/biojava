package org.biojava.nbio.structure.chem;

import org.biojava.nbio.structure.io.cif.CifBean;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Properties of a bond in a chemical component.
 * @author Sebastian Bittrich
 * @since 6.0.0
 */
public class ChemCompBond implements CifBean {
    private static final long serialVersionUID = 5905371029161975421L;
    private static final Logger logger = LoggerFactory.getLogger(ChemCompBond.class);

    private String compId;
    private String atomId1;
    private String atomId2;
    private String valueOrder;
    private String pdbxAromaticFlag;
    private String pdbxStereoConfig;
    private int pdbxOrdinal;

    public static Logger getLogger() {
        return logger;
    }

    public String getCompId() {
        return compId;
    }

    public void setCompId(String compId) {
        this.compId = compId;
    }

    public String getAtomId1() {
        return atomId1;
    }

    public void setAtomId1(String atomId1) {
        this.atomId1 = atomId1;
    }

    public String getAtomId2() {
        return atomId2;
    }

    public void setAtomId2(String atomId2) {
        this.atomId2 = atomId2;
    }

    public String getValueOrder() {
        return valueOrder;
    }

    public void setValueOrder(String valueOrder) {
        this.valueOrder = valueOrder;
    }

    public String getPdbxAromaticFlag() {
        return pdbxAromaticFlag;
    }

    public void setPdbxAromaticFlag(String pdbxAromaticFlag) {
        this.pdbxAromaticFlag = pdbxAromaticFlag;
    }

    public String getPdbxStereoConfig() {
        return pdbxStereoConfig;
    }

    public void setPdbxStereoConfig(String pdbxStereoConfig) {
        this.pdbxStereoConfig = pdbxStereoConfig;
    }

    public int getPdbxOrdinal() {
        return pdbxOrdinal;
    }

    public void setPdbxOrdinal(int pdbxOrdinal) {
        this.pdbxOrdinal = pdbxOrdinal;
    }

    /**
     * Converts this ChemCompBond's value_order attribute into an int using the
     * conversion:
     *
     * <pre>
     * 	SING -> 1
     * 	DOUB -> 2
     * 	TRIP -> 3
     * 	QUAD -> 4
     * </pre>
     *
     * Any other values will return -1.
     * <p>
     * (Source:
     * http://mmcif.rcsb.org/dictionaries/mmcif_mdb.dic/Items/_chem_comp_bond.
     * value_order.html)
     *
     * @return the numerical value of this ChemCompBond's bond order, or -1 if
     *         the value is non-numeric or unknown.
     */
    public int getNumericalBondOrder() {
        switch (valueOrder) {
            case "SING":
                return 1;
            case "DOUB":
                return 2;
            case "TRIP":
                return 3;
            case "QUAD":
                return 4;
            default:
                logger.error("Unknown or non-numeric value for value_order: " + valueOrder);
                return -1;
        }
    }
}
