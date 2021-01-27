package org.biojava.nbio.structure.chem;

import org.biojava.nbio.structure.io.cif.CifBean;

/**
 * Properties of an atom of a chemical component.
 * @author Sebastian Bittrich
 * @since 6.0.0
 */
public class ChemCompAtom implements CifBean {
    private static final long serialVersionUID = 4070599340294758941L;
    private String compId;
    private String atomId;
    private String altAtomId;
    private String typeSymbol;
    private int charge;
    private int pdbxAlign;
    private String pdbxAromaticFlag;
    private String pdbxLeavingAtomFlag;
    private String pdbxStereoConfig;
    private double modelCartnX;
    private double modelCartnY;
    private double modelCartnZ;
    private double pdbxModelCartnXIdeal;
    private double pdbxModelCartnYIdeal;
    private double pdbxModelCartnZIdeal;
    private String pdbxComponentCompId;
    private String pdbxResidueNumbering;
    private String pdbxComponentAtomId;
    private String pdbxPolymerType;
    private String pdbxRefId;
    private String pdbxComponentId;
    private int pdbxOrdinal;

    public String getCompId() {
        return compId;
    }

    public void setCompId(String compId) {
        this.compId = compId;
    }

    public String getAtomId() {
        return atomId;
    }

    public void setAtomId(String atomId) {
        this.atomId = atomId;
    }

    public String getAltAtomId() {
        return altAtomId;
    }

    public void setAltAtomId(String altAtomId) {
        this.altAtomId = altAtomId;
    }

    public String getTypeSymbol() {
        return typeSymbol;
    }

    public void setTypeSymbol(String typeSymbol) {
        this.typeSymbol = typeSymbol;
    }

    public int getCharge() {
        return charge;
    }

    public void setCharge(int charge) {
        this.charge = charge;
    }

    public int getPdbxAlign() {
        return pdbxAlign;
    }

    public void setPdbxAlign(int pdbxAlign) {
        this.pdbxAlign = pdbxAlign;
    }

    public String getPdbxAromaticFlag() {
        return pdbxAromaticFlag;
    }

    public void setPdbxAromaticFlag(String pdbxAromaticFlag) {
        this.pdbxAromaticFlag = pdbxAromaticFlag;
    }

    public String getPdbxLeavingAtomFlag() {
        return pdbxLeavingAtomFlag;
    }

    public void setPdbxLeavingAtomFlag(String pdbxLeavingAtomFlag) {
        this.pdbxLeavingAtomFlag = pdbxLeavingAtomFlag;
    }

    public String getPdbxStereoConfig() {
        return pdbxStereoConfig;
    }

    public void setPdbxStereoConfig(String pdbxStereoConfig) {
        this.pdbxStereoConfig = pdbxStereoConfig;
    }

    public double getModelCartnX() {
        return modelCartnX;
    }

    public void setModelCartnX(double modelCartnX) {
        this.modelCartnX = modelCartnX;
    }

    public double getModelCartnY() {
        return modelCartnY;
    }

    public void setModelCartnY(double modelCartnY) {
        this.modelCartnY = modelCartnY;
    }

    public double getModelCartnZ() {
        return modelCartnZ;
    }

    public void setModelCartnZ(double modelCartnZ) {
        this.modelCartnZ = modelCartnZ;
    }

    public double getPdbxModelCartnXIdeal() {
        return pdbxModelCartnXIdeal;
    }

    public void setPdbxModelCartnXIdeal(double pdbxModelCartnXIdeal) {
        this.pdbxModelCartnXIdeal = pdbxModelCartnXIdeal;
    }

    public double getPdbxModelCartnYIdeal() {
        return pdbxModelCartnYIdeal;
    }

    public void setPdbxModelCartnYIdeal(double pdbxModelCartnYIdeal) {
        this.pdbxModelCartnYIdeal = pdbxModelCartnYIdeal;
    }

    public double getPdbxModelCartnZIdeal() {
        return pdbxModelCartnZIdeal;
    }

    public void setPdbxModelCartnZIdeal(double pdbxModelCartnZIdeal) {
        this.pdbxModelCartnZIdeal = pdbxModelCartnZIdeal;
    }

    public String getPdbxComponentCompId() {
        return pdbxComponentCompId;
    }

    public void setPdbxComponentCompId(String pdbxComponentCompId) {
        this.pdbxComponentCompId = pdbxComponentCompId;
    }

    public String getPdbxResidueNumbering() {
        return pdbxResidueNumbering;
    }

    public void setPdbxResidueNumbering(String pdbxResidueNumbering) {
        this.pdbxResidueNumbering = pdbxResidueNumbering;
    }

    public String getPdbxComponentAtomId() {
        return pdbxComponentAtomId;
    }

    public void setPdbxComponentAtomId(String pdbxComponentAtomId) {
        this.pdbxComponentAtomId = pdbxComponentAtomId;
    }

    public String getPdbxPolymerType() {
        return pdbxPolymerType;
    }

    public void setPdbxPolymerType(String pdbxPolymerType) {
        this.pdbxPolymerType = pdbxPolymerType;
    }

    public String getPdbxRefId() {
        return pdbxRefId;
    }

    public void setPdbxRefId(String pdbxRefId) {
        this.pdbxRefId = pdbxRefId;
    }

    public String getPdbxComponentId() {
        return pdbxComponentId;
    }

    public void setPdbxComponentId(String pdbxComponentId) {
        this.pdbxComponentId = pdbxComponentId;
    }

    public int getPdbxOrdinal() {
        return pdbxOrdinal;
    }

    public void setPdbxOrdinal(int pdbxOrdinal) {
        this.pdbxOrdinal = pdbxOrdinal;
    }
}
